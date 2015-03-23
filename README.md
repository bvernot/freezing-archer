# freezing-archer

To run with Neanderthal as archaic:

    python ../../../bin/windowed_calculations.py \
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


To run with Denisovan as archaic:

    python ../../../bin/windowed_calculations.py \
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
