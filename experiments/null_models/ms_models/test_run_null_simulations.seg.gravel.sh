seg=$1

python ../generate_ms_params.py -mig 0 -m asn_ea_aa -N0 7310 -gen-time 25 -s2b 1032 -s3b 554 -seg $seg -recomb 1e-8 --arc-chrs 0 -p1 100 -p2 1 -p3 0 -r -ms ~/bin/ms -n10000  |  python  ../windowed_calculations.py     -vcf - -ms         -target-pops 1     -ref-pops 0     -p 10     -s-star     -winlen 50000     -winstep 50000 --ms-pop-sizes 2 200 2  --ms-num-diploid-inds 101
