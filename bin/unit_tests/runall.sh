# echo
# echo
# echo
# echo RUNNING SCRIPT
# date
# echo
# echo
# echo

baseline_git_hash=52ecafe
baseline_git_hash=32e4bc4
baseline_git_hash=189efad # local pval on S* regions
current_git_hash=$(git log --pretty=format:'%h' -n 1)

ofile_current=output_files/runall.s_star.test_$current_git_hash.txt
ofile_current_cols=output_files/runall.s_star.test_$current_git_hash.cols.txt
ofile_current_shuf=output_files/runall.s_star.test_$current_git_hash.shuf_cols.txt
ofile_current_shuf_cols=output_files/runall.s_star.test_$current_git_hash.shuf_cols.cols.txt

ofile_baseline=output_files/runall.s_star.test_$baseline_git_hash.txt
ofile_baseline_cols=output_files/runall.s_star.test_$baseline_git_hash.cols.txt
ofile_baseline_shuf=output_files/runall.s_star.test_$baseline_git_hash.shuf_cols.txt
ofile_baseline_shuf_cols=output_files/runall.s_star.test_$baseline_git_hash.shuf_cols.cols.txt

python  ../windowed_calculations.py \
    -vcf test_files/science_supplement_example.vcf \
    -indf test_files/science_supplement_example.indpop \
    -target-pops tar \
    -ref-pops ref \
    --neand-vcf test_files/science_supplement_example.arch \
    -p 10 \
    -s-star \
    -winlen 50000 \
    -winstep 50000 \
    > $ofile_current

#    | sort -k1,1 -k2,3n -k6,6 \


cat $ofile_current | transpose | awk '{print NR, $0}' > $ofile_current_cols

# awk -OFS '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $11, $12, $15, $13, $14, $10}' test_files/science_supplement_example.vcf > test_files/science_supplement_example.vcf.shuf_cols


python  ../windowed_calculations.py \
    -vcf test_files/science_supplement_example.vcf.shuf_cols \
    -indf test_files/science_supplement_example.indpop \
    -target-pops tar \
    -ref-pops ref \
    --neand-vcf test_files/science_supplement_example.arch \
    -p 10 \
    -s-star \
    -winlen 50000 \
    -winstep 50000 \
    > $ofile_current_shuf

#    | sort -k1,1 -k2,3n -k6,6 \


cat $ofile_current_shuf | transpose | awk '{print NR, $0}' > $ofile_current_shuf_cols


echo
echo "sdiff of baseline shuffled, and newest version unshuffled"
sdiff -w200 $ofile_baseline_shuf $ofile_current
echo
echo "sdiff of baseline unshuffled, and newest version shuffled (transposed)"
sdiff -w200 $ofile_baseline_cols $ofile_current_shuf_cols
