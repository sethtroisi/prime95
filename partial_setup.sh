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

    timeout "$stoptime" ./mprime -d | tee -a "$fn.log"
}


echo "Creating temporary in folder \"$1\""
mkdir "$1"
cd "$1"

ln -s ../mprime .

./mprime -v | tee version.txt

cat <<- 'EOF' > local.txt
	WorkerThreads=1
	CoresPerTest=1
EOF

cat <<- 'EOF' > prime.txt
	OutputIterations=1000000
EOF


echo -e "Y\nN\nY\n5\n" | ./mprime > /dev/null

partial_setup 4  "ECM2=1,2,14009,-1,2000000,0,5"            e0014009 "ECM, stopped at ~10% stage1"
partial_setup 5  "ECM2=1,2,14153,-1,10000,200000000,5"      e0014153 "ECM, stopped at ~20% stage2"
partial_setup 4  "ECM2=1,2,14243,-1,6000,3000000,100"       e0014243 "ECM, stopped in stage2, curve > 1"

partial_setup 5  "Pminus1=N/A,1,2,13003,-1,200000000,0"     m0013003 "PM1, Stopping in small primes ~1%"
partial_setup 25 "Pminus1=N/A,1,2,13007,-1,200000000,0"     m0013007 "PM1, Stopping in stage 1 ~7%"

partial_setup 2  "Pminus1=N/A,1,2,13009,-1,10000,10000"     m0013009 "PM1, B1 only, 1e5 finished"
partial_setup 3  "Pminus1=N/A,1,2,13121,-1,100000,90000000" m0013121 "PM1, Stopping in stage 2 ~30%"
partial_setup 2  "Pminus1=N/A,1,2,13217,-1,100000,1000000"  m0013217 "PM1, B1=1e6, B2=10e6 finished"

partial_setup 1  "PRP=1,2,500009,-1"    p0500009       "PRP, stopped at ~4%"
partial_setup 10 "PRP=46157,2,698207,1" p46157_698207  "PRP, Seventeen or Bust prime stopped at ~8%"
partial_setup 3  "PRP=6,10,71299,7"     p6_71299_7     "PRP, Near Repunit prime stopped at ~30%"

ls | grep -v '\.txt$'
cd ..
./prime95_status.py "$1"
