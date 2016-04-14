# freezing-archer

Tools for identifying introgressed archaic sequence.  These programs and scripts were used in Vernot et al, Science, 2016.  They are somewhat hacked together - please let me know if something doesn't work.

Requirements:
 - python 2.x [I use 2.7.3, so I can't guarantee that it will work with other versions]
 - python bitarray module: https://pypi.python.org/pypi/bitarray/0.8.1
 - all files from http://akeylab.gs.washington.edu/Vernot_2016/program_data

## General pipeline (details below):

1. Calculate S* and archaic match p-values in 50kb windows.
    - This is usually run once for Denisovan and once for Neandertal. Unfortunately, the code can handle only one archaic genome at a time.
2. Assign S* thresholds based on simulated data
3. Compute posterior probabilities for each putative introgressed haplotype, and categorize into null, Neandertal and Denisovan haplotypes.

## S* and Archaic match p-values

I use two custom file formats:
* .bbg - binary bed file
* .bsg - binary sequence file [essentially binary fasta]

These are my ancient attempts at having constant-time lookup to get a) if a particular site is masked or not, and b) the reference / ancestral / etc base at a given position.  These files can be generated with:

    python bin/myBedTools3.py merge -b yourbedfile.bed -obbg yourbedfile.bed.bbg

This currently only works for hg19 coordinates, and the bed file chromosomes must begin with "chr" and be in the set chr1-22,chrX,chrY,chrMT. They can be converted back with:

    python bin/myBedTools3.py merge -bbg yourbedfile.bed.bbg
    python bin/myBedTools3.py merge -bsg yourbedfile.bed.bsg

### To run with Neanderthal as archaic:

    chr=1
    arc=neand # arc=den if running for denisovan
    pop=PNG
    arc_vcf=program_data/filtered_vcfs_${arc}_mpi_minimal_filters/chr$chr.${arc}_filtered.vcf.gz
    # run on first 5mb of chromosome
    s=0
    e=5000000
    
    python  bin/windowed_calculations.py \
        --s-star \
        --vcf-has-illumina-chrnums \
        -vcfz your_data.vcf.gz \
        -indf sample_pop_mapping_file.txt \
        -ptable program_data/archaic_match_table_files.${arc}_table.fields_8-10.12-.gz.db \
        -target-pops $pop \
        -ref-pops YRI \
        --archaic-vcf $arc_vcf \
        -p 10 \
        -ancbsg program_data/chimp.bsg \
        -winlen  50000 \
        -winstep 10000 \
        -x exclude_bases.bbg \
        -r callable_bases.bbg \
        -ir intersect_callable_bases.bbg \
        -table-query mh len mapped \
        -tag-ids bigpop -tags $pop \
        -range $s $e

Format for sample_pop_mapping_file.txt:

    sample  pop     super_pop       gender
    UV043   P       PNG     .
    UV1003  AN      PNG     .
    UV1042  AN      PNG     .
    UV1043  AN      PNG     .
    UV1134  P       PNG     .


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
