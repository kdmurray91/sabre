#!/bin/bash

# exit on failure of any step
set -e

# remove output dir if it exists
test -d tmp && rm -r tmp
mkdir -p tmp
# run sabre
set -x
./sabre se -f test/test_1.fq.gz -u tmp/se_unpaired.fq -b test/bcd_se.txt
./sabre pe -f test/test_1.fq.gz -r test/test_2.fq.gz -u tmp/pe_unpaired_1.fq -w tmp/pe_unpaired_2.fq -b test/bcd_pe.txt
./sabre comb -f test/test_1.fq.gz -r test/test_c.fq.gz -u tmp/comb_unpaired_1.fq -w tmp/comb_unpaired_2.fq -b test/bcd_comb.txt
set +x
# Check output
if [ "$(find tmp/ -type f |xargs md5sum |md5sum)" == "cd6a4042a31ffc75f1c9326ed3836e0f  -" ] ; then echo "Passed!"; else echo "Failed!!!"; fi
