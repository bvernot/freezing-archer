# echo
# echo
# echo
# echo RUNNING SCRIPT
# date
# echo
# echo
# echo



python  ../windowed_calculations.py \
    -vcf science_supplement_example.vcf \
    -indf science_supplement_example.indpop \
    -target-pops tar \
    -ref-pops ref \
    --neand-vcf science_supplement_example.arch \
    -p 10 \
    -s-star \
    -winlen 50000 \
    -winstep 50000
