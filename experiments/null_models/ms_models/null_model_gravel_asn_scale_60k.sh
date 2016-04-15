seg=$1
recomb=$2
nsam=$3

module load python/2.7.3

python2 ../../bin/generate_ms_params.scale_mig.py -mig 0 -m asn_ea_aa \
    -N0 7310 -gen-time 25 -s2b 1032 \
    -s3b 554 -seg $seg -recomb ${recomb}e-8 \
    --arc-chrs 0 -p1 108 -p2 0 -p3 1 \
    -r -ms ~/bin/ms -n $nsam -split 60000 -split23 45000 \
    --rescale-migration-eur-asn-afr  \
    |  python  ../../bin/windowed_calculations.py     \
    -vcf - -ms         -target-pops 1     -ref-pops 0     \
    -p 10     -s-star     -winlen 50000     -winstep 50000 \
    --ms-pop-sizes 2 216 2  --ms-num-diploid-inds 109 \
    --tag-ids recomb --tags $recomb \
    --no-pvalues \
    --ms-simulated-region-length 50000 \
    | gzip -c 
