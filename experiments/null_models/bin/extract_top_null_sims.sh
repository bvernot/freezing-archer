for model in \
    null_simulations_null_model_two_pop_10k_5k_base \
    null_simulations_null_model_two_pop_10k_10k_base \
    ; do
    for f in $model/s_star_null.null_model_two_pop_10k_5k_base.seg_050.recomb_0.06392786121e8.iter_1.gz; model=$(basename $f | sed 's/s_star_null.//' | sed 's/\.seg.*//'); echo $model; zcat $f | sed 1d | awk 'BEGIN {OFS="\t"}  {print $4,$9,$NF,"'$model'"}' | sort -k2,2nr | head
