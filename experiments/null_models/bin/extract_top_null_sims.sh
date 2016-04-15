for model in \
    null_simulations_null_model_gravel_asn_scale_60k \
    null_simulations_null_model_gravel_eur_scale_60k \
    null_simulations_null_model_gravel_asn_scale_128k \
    null_simulations_null_model_gravel_eur_scale_128k \
    ; do

    #ofile=$model.top10k.summary
    ofile=$model.all_sims.summary.gz

    echo $model to file $ofile

    echo header should be:
    echo n_region_ind_snps       s_star  recomb

    ## $NF to get the last field (which for now is the recomb rate tag)
    header=$(zcat $model/*seg_010.recomb_0.00005829466e8.iter_1.gz \
        | awk 'BEGIN {OFS="\t"}  {print $6,$9,$NF}' \
        | head -n1)

    echo $header
    echo $header | gzip -c > $ofile

    for f in $(find $model -name '*.gz' ! -empty) ; do
        
        echo $f
        #fbase=$(basename $f .iter_1.gz)
        #echo $fbase
        #echo $model/$fbase.iter_*.gz
        #echo $f
        
        # zcat $model/$fbase.iter_*.gz \
        #     | sed 1d \
        #     | awk 'BEGIN {OFS="\t"}  {print $6,$9,$NF,"'$model'"}' \
        #     | sort -k2,2nr \
        #     | head -n10000 \
        #     >> $ofile

        ## $NF to get the last field (which for now is the recomb rate tag)
        zcat $f \
            | sed 1d \
            | awk 'BEGIN {OFS="\t"}  {print $6,$9,$NF}' \
            | gzip -c \
            >> $ofile
        
        #exit

    done

    # exit

done
