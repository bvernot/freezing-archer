# freezing-archer

Tools for identifying introgressed archaic sequence.  These programs and scripts were used in Vernot et al, Science, 2016.  They are somewhat hacked together - please let me know if something doesn't work.

## General pipeline (details below):

1. Calculate S* and archaic match p-values in 50kb windows.
    - This is usually run once for Denisovan and once for Neandertal. Unfortunately, the code can handle only one archaic genome at a time.
2. Compute posterior probabilities for each putative introgressed haplotype, and categorize into null, Neandertal and Denisovan haplotypes.

## S* and Archaic match p-values

### To run with Neanderthal as archaic:

    chr=1
    
    python bin/windowed_calculations.py \
     --vcf-has-illumina-chrnums \
     -vcfz /net/akey/vol1/home/bvernot/archaic_1kg_p3/data/png_phased_vcfs/round_3_merged_with_1kg/phased.png.chr$chr.merged_with_1kg.vcf.gz \
     -indf /net/akey/vol1/home/bvernot/archaic_1kg_p3/data/png_phased_vcfs/sample_id_files/demographics.txt.sorted.just_35_seqed_inds.language.with_1kg \
     -target-pops EUR SAS EAS PNG \
     -ref-pops YRI \
     --archaic-vcf /net/akey/vol1/home/bvernot/archaic_exome/data/neanderthal_altai_vcfs/2014.09.29/filtered_vcfs/chr$chr.altai_neand_filtered.vcf.gz \
     -p 10 \
     -s-star \
     -ancbsg /net/akey/vol1/home/bvernot/archaic_exome/data/chimp_from_den_epo/latest/chimp_chrAll.bsg \
     -winlen 50000 \
     -winstep 20000 \
     -x \
     /net/akey/vol1/home/bvernot/tishkoff/filter_files/cpg2_files/cpg2_windows_hg19.bed.zerobased.bbg \
     /net/akey/vol1/home/bvernot/archaic_exome/data/segdups/2013.01.04/genomicSuperDups.txt.bed.bbg \
     /net/akey/vol1/home/bvernot/archaic_exome/data/snp_mappability_reich/2013.01.04/hs37m_mask35_50.flt.bed.fixed.bbg \
     -r /net/akey/vol1/home/bvernot/archaic_exome/data/neanderthal_altai_vcfs/2014.09.29/neand_called_bases_x_indels.bbg \
     -ir /net/akey/vol1/home/bvernot/archaic_exome/data/chimp_from_den_epo/latest/chimp_chrAll.mapped.bbg


### To run with Denisovan as archaic:

    chr=1
    
    python bin/windowed_calculations.py \
     --vcf-has-illumina-chrnums \
     -vcfz /net/akey/vol1/home/bvernot/archaic_1kg_p3/data/png_phased_vcfs/round_3_merged_with_1kg/phased.png.chr$chr.merged_with_1kg.vcf.gz \
     -indf /net/akey/vol1/home/bvernot/archaic_1kg_p3/data/png_phased_vcfs/sample_id_files/demographics.txt.sorted.just_35_seqed_inds.language.with_1kg \
     -target-pops SAS EAS PNG \
     -ref-pop YRI \
     --archaic-vcf /net/akey/vol1/home/bvernot/archaic_exome/data/denisova_vcfs/2013.06.18/filtered_vcfs/chr$chr.den_filtered.vcf \
     -p 10 \
     -s-star \
     -ancbsg /net/akey/vol1/home/bvernot/archaic_exome/data/chimp_from_den_epo/latest/chimp_chrAll.bsg \
     -winlen 50000 \
     -winstep 20000 \
     -x \
     /net/akey/vol1/home/bvernot/tishkoff/filter_files/cpg2_files/cpg2_windows_hg19.bed.zerobased.bbg \
     /net/akey/vol1/home/bvernot/archaic_exome/data/segdups/2013.01.04/genomicSuperDups.txt.bed.bbg \
     /net/akey/vol1/home/bvernot/archaic_exome/data/snp_mappability_reich/2013.01.04/hs37m_mask35_50.flt.bed.fixed.bbg \
     -r /net/akey/vol1/home/bvernot/archaic_exome/data/denisova_vcfs/2013.06.18/denisova_called_bases_x_indels.bbg \
     -ir /net/akey/vol1/home/bvernot/archaic_exome/data/chimp_from_den_epo/latest/chimp_chrAll.mapped.bbg




### Running on simulated data [not necessary for most analyses]

#### A toy run with ms output

A toy example ms command with 4 Africans, 2 non-Africans, and one archaic chromosome:

    ms 13 1 -s 20 -I 3 8 4 1 .25

And then to calculate S* on that command:

    ms 13 1 -s 20 -I 3 8 4 1 .25 \
     | python bin/windowed_calculations.py     \
     -vcf - -ms         -target-pops 1     -ref-pops 0     \
     -p 10     -s-star     -winlen 50000     -winstep 50000 \
     --ms-pop-sizes 3 8 4 1  --ms-num-diploid-inds 6 \
     -msarc 2

A few details:
* **--ms-pop-sizes 3 8 4 1:** Formatted like ms's -I parameter
* **--ms-num-diploid-inds 6:** The number of diploid modern human individuals
* **-msarc 2:** The archaic population, numbered from 0.  This *has* to be the last population (i.e., here we're simulating populations 0,1,2)

#### A more realistic ms command

Population 1 is Africans, Population 2 is East Asian, Population 3 is European.  In bin/generate_ms_params.scale_mig.py, -p1 108 -p2 0 -p3 1 means simulate 108 Africans, no East Asians, and 1 European.  No archaic individuals are simulated.

    seg=100
    recomb=1
    nsam=10

    python2 bin/generate_ms_params.scale_mig.py -mig 0 -m asn_ea_aa \
     -N0 7310 -gen-time 25 -s2b 1032 \
     -s3b 554 -seg $seg -recomb ${recomb}e-8 \
     --arc-chrs 0 -p1 108 -p2 0 -p3 1 \
     -r -ms ~/bin/ms -n $nsam -split23 45000 \
     --rescale-migration-eur-asn-afr  \
     |  python  ../../bin/windowed_calculations.py     \
     -vcf - -ms         -target-pops 1     -ref-pops 0     \
     -p 10     -s-star     -winlen 50000     -winstep 50000 \
     --ms-pop-sizes 2 216 2  --ms-num-diploid-inds 109 \
     --tag-ids recomb --tags $recomb \
     --no-pvalues

[More models can be found here](experiments/null_models/ms_models)

[Commands to run S* on those models here](experiments/null_models/bin/submit_null_model_grid_simulations.sh)


## Assign S* Thresholds from Simulated Data

This requires a precomputed glm model, fitting recombination rate and diversity to S* quantiles.  This model is in the supporting data folder.

    for f in s_star_results*.gz ; do
        bin/process_sstar_into_haplotypes3.R $f
    done

## Compute Posterior Probabilities and Assign Introgressed Status

Combine all chromosome files into one large file for neand and one for den:
    bin/tsvcatgz s_star_results_neand_chr*.gz.processed_haps.gz | gzip -c > s_star_results_neand.ALLCHRS.gz.processed_haps.gz
    bin/tsvcatgz s_star_results_den_chr*.gz.processed_haps.gz | gzip -c > s_star_results_den.ALLCHRS.gz.processed_haps.gz

make the outputdir:
mkdir output

Run the script:

time Rscript ../bin/pval_LL_methods_pick_a_model.R RPS query_trial_output6.new_filters.tsvcat.RPS.neand.ALLCHRS.gz.processed_haps.gz query_trial_output6.new_filters.tsvcat.RPS.den.ALLCHRS.gz.processed_haps.gz output
