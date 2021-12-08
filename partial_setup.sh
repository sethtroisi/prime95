#!/usr/bin/bash

# Try to create partial files at various points
#
# Copyright 2021 Seth Troisi
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# this is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# If not, see <https://www.gnu.org/licenses/>.

set -ex

shopt -s expand_aliases

alias trace_on='set -x'
alias trace_off='{ set +x; } 2>/dev/null'

partial_setup() {
    trace_off
    stoptime="$1" # e.g. 2.5 (seconds)
    line="$2"     # e.g. "PRP=1,2,500009,-1"
    fn="$3"       # expected filename e.g. "m013003"
    comment="$4"  # e.g. "Stopping at ~7%"

    echo "$line" > "worktodo.txt"
    trace_on
    echo "$comment" > "$fn.log"

    timeout "$stoptime" ./mprime_test -d | tee -a "$fn.log"
}


if test "$#" -ne 2; then
  echo "Usage: $0 <TEST_DIR_NAME> <MPRIME_BINARY_PATH>"
  exit 1
fi

echo "Creating temporary in folder \"$2\""
mkdir "$2"
cd "$2"

echo "Using mprime from \"$1\""
if test -f "$1"; then
  ln -s "$1" mprime_test
elif test -f "../$1"; then
  ln -s "../$1" mprime_test
else
  echo "can't find mprime at \"$1\""
  exit 1
fi

./mprime_test -v | tee version.txt

cat <<- 'EOF' > local.txt
	WorkerThreads=1
	CoresPerTest=1
EOF

cat <<- 'EOF' > prime.txt
	OutputIterations=5000000
EOF

# Join Gimps, Use Primenet, Accept, Settings (default, default, default), exit
echo -e "Y\nN\nY\n\n\n\n\n\n\n\n\n\n\n\n\n5\n" | ./mprime_test > /dev/null

partial_setup 4  "ECM2=1,2,14009,-1,2000000,0,5"            e0014009 "ECM, stopped at ~10% stage1"
# Hopefully stops in stage3 but can't guarentee
partial_setup 4  "ECM2=1,2,14243,-1,50000,3000000,100"      e0014243 "ECM, likely in stage2, curve > 1"
# Large M makes it easier to stop in stage 2
partial_setup 19 "ECM2=1,2,150089,-1,50000,30000000,5"      e0150089 "ECM, stopped at ~40% stage2"

partial_setup 8  "Pminus1=N/A,1,2,2237,-1,200000000,0"      m0002237 "PM1, Stopping in small primes ~2%"
# Take us to close to the threshold
partial_setup 100 "Pminus1=N/A,1,2,2267,-1,3000000000,0"      m0002267 "PM1, Stopping in stage 1 ~20%"
# Take us over the threshold (with a new save file)
partial_setup 40 "Pminus1=N/A,1,2,2267,-1,3000000000,0"      m0002267 "PM1, Stopping in stage 1 ~20%"

partial_setup 2  "Pminus1=N/A,1,2,13009,-1,100000,100000"     m0013009 "PM1, B1 only, 1e5 finished"
partial_setup 2  "Pminus1=N/A,1,2,13217,-1,100000,1000000"  m0013217 "PM1, B1=1e5, B2=1e6 finished"
# Large M makes it easier to stop in stage 2
partial_setup 20  "Pminus1=N/A,1,2,150107,-1,100000,4000000000" m0150107 "PM1, Stopping in stage 2 ~20%"

partial_setup 2  "PRP=1,2,700001,-1"    p0700001       "PRP, stopped at ~3%"
partial_setup 10 "PRP=46157,2,698207,1" p46157_698207  "PRP, Seventeen or Bust prime stopped at ~8%"
partial_setup 3  "PRP=6,10,71299,7"     p6_71299_7     "PRP, Near Repunit prime stopped at ~30%"

ls | grep -v '\.txt$'
cd ..
./prime95_status.py "$2" --json "$2.json"
