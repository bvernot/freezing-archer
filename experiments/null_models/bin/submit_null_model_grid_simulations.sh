echo
echo
echo "=================================================="
echo "           SUBMITTING NULL GRID SIMS"
echo "=================================================="
echo
echo


# recombination rates to simulate (from science paper)
recomb_grid="0.00005829466
  0.00009611165
  0.00015846133
  0.00026125856
  0.00043074254
  0.00071017439
  0.00117087962
  0.00193045414
  0.00318278080
  0.00524751840
  0.00865169520
  0.01426423391
  0.02351774586
  0.03877420783
  0.06392786121
  0.10539922456
  0.17377394345
  0.28650479686
  0.47236655274
  0.77880078307
  1.28402541669
  2.11700001661
  3.49034295746
  5.75460267601
  9.48773583636
 15.64263188419"

recomb_grid="0.00005829466
  0.00009611165
  0.00015846133
  0.00026125856
  0.00043074254
  0.00071017439
  0.00117087962
  0.00193045414
  0.00318278080
  0.00524751840
  0.00865169520
  0.01426423391
  0.02351774586
  0.03877420783
  0.06392786121
  0.10539922456
  0.17377394345
  0.28650479686
  0.47236655274
  0.77880078307
  1.28402541669
  2.11700001661
  3.49034295746
  5.75460267601
  9.48773583636
 15.64263188419"

# for recomb in $recomb_grid ; do
#     echo "|${recomb}|"
# done



for model in ms_models/null_model_gravel_asn_scale_60k.sh \
    ms_models/null_model_two_pop_10k_10k_base.sh \
    ms_models/null_model_two_pop_10k_5k_base.sh ; do
    
    iter=1
    n=50
    n=50000
    model_name=$(basename $model .sh)
    
    if [ ! -e null_simulations_$model_name ] ; then
        echo "making null directory: null_simulations_$model_name"
        mkdir null_simulations_$model_name
    fi
    
    for seg in $(seq -f "%03g" 10 5 400) ; do

        for recomb in $(echo $recomb_grid | cut -f 1-3 -d ' ') ; do
        #for recomb in 0.77880078307 ; do

            ofile=null_simulations_$model_name/s_star_null.$model_name.seg_$seg.recomb_${recomb}e8.iter_$iter.gz

            if [ ! -e $ofile ] ; then
                
                echo $iter $seg $recomb $model_name

                touch $ofile
                
                qsub -l mfree=.1G -js 0 /net/akey/vol1/home/bvernot/tishkoff/run_cluster_job.save_output.sh `pwd` null $ofile \
                    bash $model $seg $recomb $n

            fi
            
            #break
            
        done
        
        #break

    done

done


exit

-------------------------
 Sun Mar 15 14:33:22 PDT 2015 
-------------------------

time for f in test_null_simulations/test.s_star_null.seg_*.recomb_1e8.iter_1.txt ; do sed 1d $f | cut -f 4,5,8,9 ; done > test_null_simulations.summary.seg

time for f in test_null_simulations/test.s_star_null.mu2_* ; do mu=$(echo $f | sed 's/.*mu2_//g' | sed 's/.recomb.*//g' | sed 's/^0//g'
); sed 1d $f | cut -f 4,5,8,9 | awk 'BEGIN{OFS="\t"} {print $0, '$mu'}' ; done > test_null_simulations.summary.mu

^^ this fixed a problem where awk interpreted 025 as 21, FOR SOME INCREDIBLY FRUSTRATING REASON.


-------------------------
 Sun Mar 15 14:45:16 PDT 2015 
-------------------------

[a002 bin/unit_tests]$ for i in {1..5}; do for seg in $(seq -f "%03g" 25 5 400) ; do echo $i $seg seg.gravel; qsub -l mfree=.05G  /net/akey/vol1/home/bvernot/tishkoff/run_cluster_job.save_output.sh `pwd` null test_null_simulations/test.s_star_null.seg.gravel_$seg.recomb_1e8.iter_$i.txt bash test_run_null_simulations.seg.gravel.sh $seg; done; done >> jobs.txt                                                                                                                 



-------------------------
 Sun Mar 15 15:04:23 PDT 2015 
-------------------------

[a002 bin/unit_tests]$ for i in {1..5}; do for seg in $(seq -f "%03g" 25 5 400) ; do echo $i $seg seg.gravel; qsub -l mfree=.05G  /net/akey/vol1/home/bvernot/
tishkoff/run_cluster_job.save_output.sh `pwd` null test_null_simulations/test.s_star_null.seg.gravel_$seg.recomb_1e8.iter_$i.txt bash test_run_null_simulation
s.seg.gravel.sh $seg; done; done >> jobs.txt                                                                                                                  
[a002 bin/unit_tests]$ for i in {1..5}; do for seg in $(seq -f "%03g" 25 5 400) ; do echo $i $seg seg.gravel.asn; qsub -l mfree=.05G  /net/akey/vol1/home/bver
not/tishkoff/run_cluster_job.save_output.sh `pwd` null test_null_simulations/test.s_star_null.seg.gravel.asn_$seg.recomb_1e8.iter_$i.txt bash test_run_null_si
mulations.seg.gravel.asn.sh $seg; done; done >> jobs.txt                                                                                                      


-------------------------
 Sun Mar 15 17:22:39 PDT 2015 
-------------------------

[a002 bin/unit_tests]$ time for f in test_null_simulations/test.s_star_null.seg*.recomb_1e8.iter_*.txt ; do cat=$( echo $f | sed 's/.*null.//g' | sed 's/_.*//g' ); sed 1d $f | cut -f 4,5,8,9 | awk 'BEGIN{OFS="\t"} {print $0, "'$cat'"}' ; done > test_null_simulations.summary.seg.cats


