seg=$1
recomb=$2
nsam=$3

module load python/2.7.3

python2 ../../bin/generate_ms_params.scale_mig.py -m two_pop -mig 0 \
    -arc-chrs 0 -s1 10000 -s2 5000 -s12 10000 -split 60000 \
    -gen-time 25 -seg $seg -recomb ${recomb}e-8 --arc-chrs 0 \
    -p1 108 -p2 1 -r -ms ~/bin/ms -n$nsam \
    | python  ../../bin/windowed_calculations.py     -vcf - -ms         \
    -target-pops 1     -ref-pops 0     -p 10     \
    -s-star     -winlen 50000     -winstep 50000 \
    --ms-pop-sizes 2 216 2  --ms-num-diploid-inds 109 \
    --tag-ids recomb --tags $recomb \
    --no-pvalues \
    | gzip -c 
