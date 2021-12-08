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

"""
This has been tested with
  * v29.8 build 6
  * v30.3 build 6
  * v30.7 build 9
  * v30.8 build 4 (partial)
"""


import argparse
import json
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
 # Removed {ECM,PRP,PM1}_VERSION
FACTOR_VERSION          = 1
LL_VERSION              = 1
SPOOL_FILE_VERSION      = 1

##### END MAGIC NUMBERS             #####


MAX_BACKUP_FILES   = 100
BACKUP_CWD         = "Status of files in '%s'."
BACKUP_CWD_ERROR   = "Unable to read working directory."
BACKUP_STATUS      = "Backup %-16s | %s."
BACKUP_NONE        = "No Backup files (*.bu) were found in %s."
BACKUP_PARSE_ERROR = "Unable to parse (%s)."

BACKUP_PTN = re.compile("[emp][0-9]+(_[0-9]+){0,2}(.bu[0-9]*)?$")



def get_arg_parser():
    parser = argparse.ArgumentParser(description="Parse prime95 status/backup files")

    parser.add_argument('dir', type=str, default=".",
        help="Directory for with status/backup files")

    parser.add_argument('--json', type=str, default="",
        help="Save JSON data in this file")

    parser.add_argument('--skip-failed', action="store_true",
        help="Don't output anything for failed files")

    return parser

def scan_directory(dir_name):
    names = []
    if not os.path.isdir(dir_name):
        sys.exit(f"{dir_name!r} does not exist")

    for filename in os.listdir(dir_name):
        if BACKUP_PTN.match(filename):
            names.append(filename)

    return names


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

def read_uint64(f):
    """Read a uint64_t from [f]ile"""
    return _read_struct(f, 8, "Q")

def read_double(f):
    """Read a double from [f]ile"""
    return _read_struct(f, 8, "d")

def read_int(f):
    """Read a slong then convert to int[32]"""
    return read_slong(f) & 0xFFFFFFFF

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

    wu["checksum"] = read_long(f)

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
            wu["work_type"] = "WORK_ECM"

            if version <= 2:    # 25 - 30.6
                wu["state"] = read_long(f)
                wu["curve"] = read_long(f)    # 'curves_to_go' in older code
                wu["sigma"] = read_double(f)  # 'curve' in older code

                wu["B"] = read_uint64(f)
                wu["B_done"] = read_uint64(f)
                wu["C_done"] = read_uint64(f)

                wu["stage_guess"] = wu["state"] + 1

            elif version == 3:
                wu["curve"] = read_long(f)
                wu["average_B2"] = read_uint64(f)
                state = read_int(f)
                wu["state"] = state

                wu["sigma"] = read_double(f)
                wu["B"] = read_uint64(f)
                wu["C"] = read_uint64(f)

                #define ECM_STATE_STAGE1_INIT       0   /* Selecting sigma for curve */
                #define ECM_STATE_STAGE1        1   /* In middle of stage 1 */
                #define ECM_STATE_MIDSTAGE      2   /* Stage 2 initialization for the first time */
                #define ECM_STATE_STAGE2        3   /* In middle of stage 2 (processing a pairmap) */
                #define ECM_STATE_GCD           4   /* Stage 2 GCD */

                if state == 1:    # ECM_STATE_STAGE1
                    wu["stage_guess"] = 1
                    wu["stage1_prime"] = read_uint64(f)
                elif state == 2:  # ECM_STATE_MIDSTAGE
                    wu["stage_guess"] = 1
                elif state == 3:  # ECM_STATE_STAGE2
                    wu["stage_guess"] = 2
                    # A bunch of unused (by this program values)
                    # read 6 ints (stage2_numvals, totrels, D, E, two_fft_stage2), pool_type)
                    # read 2 uint64 (first_relocatable, last_relocatable)
                    for i in range(6):
                        read_int(f)
                    for i in range(2):
                        read_uint64(f)
                    wu["B2_start"] = read_uint64(f)
                    wu["C_done"] = read_uint64(f)
                    pct = (wu["C_done"] - wu["B2_start"]) / (wu["C"] - wu["B2_start"])
                    wu["pct_guess"] = f"~~{pct:.1%}"
                elif state == 4:  # ECM_STATE_CGD
                    wu["stage_guess"] = 2

            else:
                sys.stderr.write(f"ECM with version {version} {filename!r}\n")
                return None

        elif magic == PM1_MAGICNUM:
            wu["work_type"] = "WORK_PMINUS1"

            if 5 <= version <= 7:    # 30.4 to 30.7
                #define PM1_STATE_STAGE0	0	/* In stage 1, computing 3^exp using a precomputed mpz exp */
                #define PM1_STATE_STAGE1	1	/* In stage 1, processing larger primes */
                #define PM1_STATE_MIDSTAGE	2	/* Between stage 1 and stage 2 */
                #define PM1_STATE_STAGE2	3	/* In middle of stage 2 (processing a pairmap) */
                #define PM1_STATE_GCD		4	/* Stage 2 GCD */
                #define PM1_STATE_DONE		5	/* P-1 job complete */

                state = read_int(f)
                wu["state"] = state

                if state == 0:    # PM1_STATE_STAGE0
                    wu["stage_guess"] = "B1_pre"
                    wu["interim_B"] = read_uint64(f)
                    wu["max_stage0_prime"] = read_long(f)
                    wu["stage0_bitnum"] = read_long(f)

                    if version == 7:
                        # TODO verify this after 30.8 stable
                        wu["stage0_bitnum"] = read_long(f)

                    wu["B1_guess"] = wu["stage0_bitnum"]

                elif state == 1:  # PM1_STATE_STAGE1
                    wu["stage_guess"] = "B1"
                    wu["B_done"] = read_uint64(f)
                    wu["interim_B"] = read_uint64(f)
                    wu["stage1_prime"] = read_uint64(f)

                    # stage1_prime can be slightly larger than bound
                    # use B_done if non-zero otherwise stage1_prime
                    wu["B1_guess"] = wu["B_done"] or wu["stage1_prime"]

                elif state == 2:  # PM1_STATE_MIDSTAGE
                    wu["stage_guess"] = "B1"
                    wu["B_done"] = read_uint64(f)
                    wu["C_done"] = read_uint64(f)

                    # TODO can C_done be > B_done if B2 was extended?
                    wu["B1_guess"] = wu["B_done"]

                elif state == 3:  # PM1_STATE_STAGE2
                    wu["stage_guess"] = "B2"
                    wu["B_done"] = read_uint64(f)
                    wu["C_done"] = read_uint64(f)
                    wu["interim_C"] = read_uint64(f)
                    wu["pct_guess"] = wu["interim_C"] / wu["C_done"]

                    wu["B1_guess"] = wu["B_done"]
                    wu["B2_guess"] = wu["C_done"]

                elif state == 4:  # PM1_STATE_GCD
                    wu["stage_guess"] = "B2"
                    wu["pct_guess"] = 0.99
                    wu["B_done"] = read_uint64(f)
                    wu["C_done"] = read_uint64(f)

                    wu["B1_guess"] = wu["B_done"]
                    wu["B2_guess"] = wu["C_done"]

                elif state == 5:  # PM1_STATE_DONE
                    wu["stage_guess"] = "DONE"
                    wu["pct_complete"] = 1
                    wu["B_done"] = read_uint64(f)
                    wu["C_done"] = read_uint64(f)

                    wu["B1_guess"] = wu["B_done"]
                    wu["B2_guess"] = wu["C_done"]

            elif version < 5:  # Version 25 through 30.3 save file
                state = read_long(f)
                wu["state"] = state
                # PM1_STATE enum changed so we instead store stage below as B1 or B2

                if version <= 2:
                    wu["max_stage0_prime"] = 13333333
                else:
                    wu["max_stage0_prime"] = read_long(f)

                # /* Read the first part of the save file, much will be ignored
                #    but must be read for backward compatibility */

                wu["B_done"]  = read_uint64(f)
                wu["B"]       = read_uint64(f)
                wu["C_done"]  = read_uint64(f)

                wu["C_start_unused"] = read_uint64(f)
                wu["C_unused"]       = read_uint64(f)  # C_done in source code, but I think this is C actually

                # "Processed" is number of bits in state 0, number of primes in state 1
                wu["processed"] = read_uint64(f)

                wu["D"]       = read_long(f)
                wu["E"]       = read_long(f)
                wu["rels_done"] = read_long(f)

                # /* Depending on the state, some of the values read above are not meaningful. */
                # /* In stage 0, only B and processed (bit number) are meaningful. */
                # /* In stage 1, only B_done, B, and processed (prime) are meaningful. */
                # /* In stage 2, only B_done is useful.  We cannot continue an old stage 2. */
                # /* When done, only B_done and C_done are meaningful. */
                if state == 3:     # PM1_STATE_STAGE0
                    wu["stage_guess"] = "B1_pre"
                    wu["B1_guess"] = wu["processed"]
                    if version == 1:
                        # 29.4 build 7 changed the calc_exp algorithm and invalidated this
                        wu["B1_guess"] = 0
                elif state == 0:  # PM1_STATE_STAGE1
                    wu["stage_guess"] = "B1"
                    wu["B1_guess"] = max(wu["B_done"], wu["processed"])
                elif state == 1:  # PM1_STATE_STAGE2
                    # Cannot continue stage 2 from old P-1 save file (so through away that data)
                    wu["stage_guess"] = "B2"
                    wu["B1_guess"] = wu["B_done"]
                    wu["B2_guess"] = wu["C_done"]
                elif state == 2:  # PM1_STATE_DONE
                    wu["stage_guess"] = "DONE"
                    wu["B1_guess"] = wu["B_done"]
                    wu["B2_guess"] = wu["C_done"]

            else:
                sys.stderr.write(f"P-1 with version {version} {filename!r}\n")
                return None


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
            if version > 7:
                sys.stderr.write(f"PRP with version {version} {filename!r}\n")
                return None

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
            sys.stderr.write(f"Unknown type magicnum = {magic}\n")
            return None

    return wu


def one_line_status(fn, wu, name_pad):
    buf = ""
    # TODO should I print {k} * {b} ^ {n} - {c} ?

    work = wu["work_type"]
    if work == "WORK_ECM":
        # TODO print bounds?
        pct = wu.get("pct_complete", wu.get("pct_guess", ""))
        if isinstance(pct, (int, float)):
            pct = f"{pct:.1%}"
        buf += "ECM | Curve {:d} | Stage {}".format(wu["curve"], wu["stage_guess"])
        if pct:
            buf += " (" + pct + ")"
    elif work == "WORK_PMINUS1":
        stage = wu["stage_guess"]
        if stage == "B1_pre":
            # Stage 1, processed = bit_number
            buf += "P-1 | Stage 1 ({:.1%}) B1 <= {:d}".format(
                wu["pct_complete"], wu["B1_guess"])
        elif stage == "B1":
            # Stage 1 after small primes
            if wu["pct_complete"] == 1:
                buf += "P-1 | B1={:d} complete".format(wu["B1_guess"])
            else:
                buf += "P-1 | Stage 1 ({:.1%}) B1 @ {:d}".format(
                    wu["pct_complete"], wu["B1_guess"])
        elif stage == "B2":
            # Stage 2
            buf += "P-1 | B1={:.0f} complete, Stage 2".format(wu["B2_guess"])
            if wu["pct_complete"] == 0.99:
                buf += " (99%, computing GCD)"
            else:
                buf += " ({:.1%})".format(wu["pct_complete"])
        elif stage == "DONE":
            # P-1 done
            buf += "P-1 | B1={:.0f}".format(wu["B1_guess"])
            if wu["B2_guess"] > wu["B1_guess"]:
                buf += ", B2={:.0f}".format(wu["B2_guess"])
                if wu.get("E", 0) >= 2:
                    buf += ", E={:d}".format(wu["E"])
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

    return fn.ljust(name_pad) + " | " + buf


def main(args):
    names = sorted(scan_directory(args.dir))

    parsed = {}
    failed = []
    for name in names:
        result = parse_work_unit_from_file(os.path.join(args.dir, name))
        if result is not None:
            parsed[name] = result
        else:
            failed.append(name)

    if failed and not args.skip_failed:
        print()
        print(f"FAILED ({len(failed)}):")
        for i, name in enumerate(failed):
            if i < 10 or (i < 100 and i % 10 == 0) or i % 100 == 0: print(f"\t{i} {name}")
        print()

    longest_name = min(20, max(map(len, parsed.keys()), default=0))
    print(f"Found {len(names)} backup files in {args.dir!r}")
    for name in sorted(parsed):
        print(one_line_status(name, parsed[name], longest_name))

    if args.json:
        print(f"Writing json data to {args.json!r}")

        # stage clutters json so pop it
        for wu in parsed.values():
            wu.pop("stage", None)

        with open(args.json, "w") as f:
            json.dump(parsed, f, indent=4)



if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()

    main(args)
