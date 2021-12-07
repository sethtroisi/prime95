#!/usr/bin/env python3

# Copyright (c) 2021 Seth Troisi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
Read Prime95 backup files and display status on them

An attempt to merge this into Prime95's c code was made but there was little
interest and the code was verbose and difficult to maintain
See
  * https://www.mersenneforum.org/showthread.php?p=540022#post540022
  * https://github.com/sethtroisi/prime95/pull/1


This was written looking at 29.8 source code but the save file header are
likely to be fairly consistent and can be upgraded easily.
The relevant files for details are
  * commonc.c
    read_header / write_header
    read_X

  * ecm.c
    pm1_save
    calc_exp
"""

"""
ostensible the file format is
    u32             magic number  (different for ll, p-1, prp, tf, ecm)
    u32             version number
    double          k in k*b^n+c
    u32             b in k*b^n+c
    u32             n in k*b^n+c
    s32             c in k*b^n+c
    double          pct complete
    char(11)        stage
    char(1)         pad
    u32             checksum of all following data
"""


import os
import re
import struct
import sys


##### MAGIC NUMBERS USED BY PRIME95  #####

 #### $ grep '#define.*MAGICNUM' *.c #####
FACTOR_MAGICNUM         = 0x1567234D
LL_MAGICNUM             = 0x2c7330a8
PRP_MAGICNUM            = 0x87f2a91b
SPOOL_FILE_MAGICNUM     = 0x73d392ac
ECM_MAGICNUM            = 0x1725bcd9
PM1_MAGICNUM            = 0x317a394b
 #### $ grep '#define.*VERSION' *.c #####
FACTOR_VERSION          = 1
LL_VERSION              = 1
PRP_VERSION             = 4
SPOOL_FILE_VERSION      = 1
ECM_VERSION     = 1
PM1_VERSION     = 2
 #### $ grep '#define.*_STAGE' *.c  #####
ECM_STAGE1      = 0
ECM_STAGE2      = 1

PM1_STAGE0      = 3       # In stage 1, squaring small primes
PM1_STAGE1      = 0       # In stage 1, processing larger primes
PM1_STAGE2      = 1       # In stage 2

# MANUALLY ADDED
PM1_DONE        = 2

##### END MAGIC NUMBERS             #####


MAX_BACKUP_FILES   = 100
BACKUP_CWD         = "Status of files in '%s'."
BACKUP_CWD_ERROR   = "Unable to read working directory."
BACKUP_STATUS      = "Backup %-16s | %s."
BACKUP_NONE        = "No Backup files (*.bu) were found in %s."
BACKUP_PARSE_ERROR = "Unable to parse (%s)."

BACKUP_PTN = re.compile("[emp][0-9]+(.bu[0-9]*)?$")




def scan_directory(dir_name):
    names = []
    if not os.path.isdir(dir_name):
        sys.exit(f"{dir_name!r} does not exist")

    for filename in os.listdir(dir_name):
        if BACKUP_PTN.match(filename):
            names.append(filename)

    return names


# Global "sum" used as checksum
global_checksum = 0

def _read_bytes(f, count):
    f_bytes = f.read(count)
    if len(f_bytes) != count:
        sys.stderr.write(f"Not enough bytes! Read {len(f_bytes)} of {count} from {f.name!r}")
        return b"\0"
    return f_bytes

def _read_struct(f, count, struct_format):
    f_bytes = _read_bytes(f, count)
    tmp = struct.unpack(struct_format, f_bytes)
    assert len(tmp) == 1, tmp

    # TODO checksumming
    return tmp[0]


def read_long(f):
    """Read a uint32_t from [f]ile"""
    return _read_struct(f, 4, "I")

def read_slong(f):
    """Read a int32_t from [f]ile"""
    return _read_struct(f, 4, "i")

def read_longlong(f):
    """Read a uint64_t from [f]ile"""
    return _read_struct(f, 8, "Q")

def read_double(f):
    """Read a double from [f]ile"""
    return _read_struct(f, 8, "d")

def read_array(f, length):
    b = _read_bytes(f, length)
    # TODO checksumming
    return b

def read_header(f, wu):
    # read_header normally validates k,b,n,c which we don't do

    wu["version"] = read_long(f)
    wu["k"] = read_double(f)
    wu["b"] = read_long(f)
    wu["n"] = read_long(f)
    wu["c"] = read_slong(f)
    wu["stage"] = list(read_array(f, 11))
    wu["pad"] = list(read_array(f, 1))
    wu["pct_complete"] = read_double(f)

    wu["sum"] = read_long(f)

    wu["stage"][10] = 0
    wu["pct_complete"] = max(0, min(1, wu["pct_complete"]))

    return True


def parse_work_unit_from_file(filename):
    wu = {}

    wu["work_type"] = None

    with open(filename, "rb") as f:
        wu["magicnumber"] = magic = read_long(f)

        # Common file header
        read_header(f, wu)
        version = wu["version"]

        if magic == ECM_MAGICNUM:
            if version != ECM_VERSION:
                sys.exit(f"ECM({magic}) with version {version}!")

            wu["work_type"] = "WORK_ECM"
            wu["pm1_stage"] = read_long(f)
            wu["curves_to_do"] = read_long(f)
            # Sigma of current curve
            wu["curve"] = read_double(f)

            wu["B"] = read_longlong(f)
            wu["B_done"] = read_longlong(f)
            wu["C_done"] = read_longlong(f)

        elif magic == PM1_MAGICNUM:
            if version in (3, 6):
                # new 30.8 verson
                wu["work_type"] = "WORK_PMINUS1"

            elif version != PM1_VERSION:
                sys.exit(f"P-1({magic}) with version {version}!")

            else:
                wu["work_type"] = "WORK_PMINUS1"
                wu["pm1_stage"] = read_long(f)
                wu["B_done"]  = read_longlong(f)
                wu["B"]       = read_longlong(f)
                wu["C_done"]  = read_longlong(f)
                wu["C_start"] = read_longlong(f)
                wu["C"]       = read_longlong(f)

                # "Processed" is number of bits in stage 0, number of primes in stage 1
                wu["pairs_done"] = read_longlong(f)
                wu["D"]       = read_long(f)
                wu["E"]       = read_long(f)

        elif magic == LL_MAGICNUM:
            if version != LL_VERSION:
                sys.exit(f"LL({magic}) with version {version}!")

            # duplicated from commonb.c readLLSaveFile (minus reading data)
            wu["work_type"] = "WORK_TEST"

            # error_count is stored in E,
            # count (iterations) is stored in C,
            wu["E"] = read_long(f)
            wu["C"] = read_long(f)

        elif magic == PRP_MAGICNUM:
            if version != PRP_VERSION:
                sys.exit(f"PRP({magic}) with version {version}!")

            wu["work_type"] = "WORK_PRP"

            # error_count is stored in E,
            # count (iterations) is stored in C,
            wu["E"] = read_long(f)
            wu["C"] = read_long(f)

        elif magic == FACTOR_MAGICNUM:
            if version != FACTOR_MAGICNUM:
                sys.exit(f"FACTOR({magic}) with version {version}!")

            wu["work_type"] = "WORK_FACTOR"
            # TODO: implement WORK_FACTOR report

        else:
            sys.stderr(f"Unknown type magicnum = {magic}")
            return None

    return wu


def one_line_status(fn, wu):
    buf = ""
    # TODO should I print {k} * {b} ^ {n} - {c} ?

    work = wu["work_type"]
    if work == "WORK_ECM":
        # TODO print bounds?
        buf += "ECM | Curve {:d} | Stage {} ({:.1%})".format(
                wu["curves_to_do"], wu["pm1_stage"] + 1, wu["pct_complete"])
    elif work == "WORK_PMINUS1":
        if wu["version"] in (3, 6):
            # new 30.8 verson
            buf += " P-1 | *unhandled 30.8 data* waiting for beta "
        else:
            stage = wu["pm1_stage"]
            if stage == PM1_STAGE0:
                # Stage 1, pairs_done = processed = bit_number
                buf += "P-1 | Stage 1 ({:.1%}) B1 <{:d}".format(
                    wu["pct_complete"], wu["pairs_done"])
            elif stage == PM1_STAGE1:
                # Stage 1 after small primes, pairs_done = processed = prime
                buf += "P-1 | Stage 1 ({:.1%}) B1 @ {:d}".format(
                    wu["pct_complete"], wu["pairs_done"])
            elif stage == PM1_STAGE2:
                # Stage 2 after small primes, pairs_done = processed = B1 bound
                buf += "P-1 | B1={:.0f} complete, Stage 2 ({:.1%})".format(
                       wu["B"], wu["pct_complete"]);
            elif stage == PM1_DONE:
                # P-1 done
                buf += "P-1 | B1={:.0f}".format(wu["B"])
                if wu["C"] > wu["B"]:
                    buf += " ,B2={:.0f}".format(wu["C"])
                    if wu["E"] >= 2:
                        buf += " ,E={:d}".format(wu["E"])
                buf += " complete"
            else:
                buf += "UNKNOWN STAGE={:d}".format(stage)

    elif work == "WORK_TEST":
        buf += "LL  | Iteration {}/{} [{:0.2%}]".format(
            wu["C"], wu["n"], wu["pct_complete"])

    elif work == "WORK_PRP":
        buf += "PRP | Iteration {}/{} [{:0.2%}]".format(
            wu["C"], wu["n"], wu["pct_complete"])

    elif work == "WORK_FACTOR":
        buf += "FACTOR | *unhandled*"

    else:
        buf += "UNKNOWN work={}".format(work)

    return fn.ljust(12) + " | " + buf


def main(dir_name = "."):
    names = scan_directory(dir_name)

    parsed = {}
    failed = []
    for name in names:
        result = parse_work_unit_from_file(os.path.join(dir_name, name))
        if result is not None:
            parsed[name] = result
        else:
            failed.append(name)

    print(f"Found {len(names)} backup files in {dir_name!r}")
    for name in sorted(parsed):
        print(one_line_status(name, parsed[name]))

    if failed:
        print()
        print("FAILED:")
        for name in sorted(failed):
            print("\t" + name)
        print()

if __name__ == "__main__":
    directory = "." if len(sys.argv) < 2 else sys.argv[1]
    main(directory)
