library(data.table)

if (length(commandArgs(T)) != 1) {
   cat("Requires a gzipped archaic S* output file.\n")
   quit()
}

neand_file = commandArgs(T)[1]
cat("Processing Neanderthal file:", neand_file, "\n")

use_precomputed_model = T

if (use_precomputed_model) {
  ## loads the object G2s
  cat('Loading precomputed null model...\n')
  load('program_data/gravel_asn_scale_60k.model.Rdata')

} else {

  cat('Loading null simulations...\n')
  null.sample = fread('gunzip -c /net/akey/vol2/bvernot/archaic_1kg_p3/experiments/null_models/null_simulations_null_model_FOUR_MODEL_QUANTILES.txt.gz')
  
  cat('Building GLM model...\n')
  G2s <- gam(s_star ~ te(snps, lr, q, k=8), data = null.sample[tag == 'gravel_asn_scale_60k'], method="GCV.Cp")
  null.sample[, p.s2.gv := predict(G2s, null.sample)]
}

cat('Reading S* files...\n')
line_limit = -1
png.neand = fread(sprintf('gunzip -c %s', neand_file), nrows=line_limit)
setkey(png.neand, chrom, winstart, winend, ind_id)

# to remove duplicates..  which we have to do because I foolishly have my 
## python code compute stats for windows at both the beginning and end of 
## the specified range.  also, at the beginning of the range I don't read
## snps *before the range*, so for the window 960000-1010000 and a range 
## that starts at 1m, only snps after 1m are counted.  this is not a problem
## at the end of the range.  so here I'm removing the window with fewer
## snps, then running unique to get the remaining dups (where all the snps
## in a window are after the beginning of the range anyway).
## 
cat('Removing duplicates due to overlapping windows...\n')
nrow(unique(png.neand))
png.neand[, max_snps_for_dups := max(n_snps), list(chrom, winstart, winend, ind_id)]
png.neand = png.neand[n_snps == max_snps_for_dups]
png.neand = unique(png.neand)
nrow(png.neand)


cat('Loading and merging recombination rates...\n')
# recomb = fread('/net/akey/vol1/home/bvernot/archaic_exome/data/recombination_rates/genetic_map_HapMapII_GRCh37/2013.01.26/windows.50000.10000.bed.recomb')
recomb = fread('program_data/genetic_map_HapMapII_GRCh37/2013.01.26/windows.50000.10000.bed.recomb')
setnames(recomb, c('chrom', 'winstart', 'winend', 'covered_segments', 'covered_bases', 'recomb', 'flag'))
recomb[, chrom := suppressWarnings(as.numeric(substr(chrom, 4,6)))]
setkey(recomb, chrom, winstart, winend)
setkey(png.neand, chrom, winstart, winend)
png.neand = recomb[png.neand]


cat('Filtering regions that do not pass our filters...\n')
## filter png.neand0
png.neand[, filter := (flag != '.' | recomb == 0 | callable_bases < 25000 | n_region_ind_snps < 100)]
write.table(png.neand[filter != T, .N, keyby=list(chrom, winstart, winend)], sprintf('%s.callable_windows.bed', neand_file), sep='\t', quote=F, col.names=T, row.names=F)

png.neand[, filter := (filter | s_star <= 0)]
png.neand = png.neand[filter != T]



ss_quantiles_fixed3 = c(.9, .95, .96, .97, .98, .985, .99, .9925, .995, .9975, .999, .9995)
ss_quantiles_fixed3 = c(.99)
ss_quantiles_fixed3 = c(.90, .91, .92, .93, .94, .95, .96, .97, .98, .985, .99, .995, .999, .9999)
cat("Calculating S* significant values...\n")
for ( qt in ss_quantiles_fixed3 ) {
  cat(qt, "\n")
  png.neand[, paste0('p.gv', qt) := predict(G2s, png.neand[, list(snps = n_region_ind_snps, lr = log(recomb), q=qt)])]
  qt
}
ss_quantiles_fixed3 = sort(ss_quantiles_fixed3)
ss_colnames_fixed3 = paste0('p.gv', ss_quantiles_fixed3)



cat('Setting the quantile for each row...\n')
png.neand[, ss_sig := 0.0]
for (i in 1:length(ss_colnames_fixed3)) {
  png.neand[eval(parse(text=sprintf("s_star >= %s", ss_colnames_fixed3[i]))), ss_sig := ss_quantiles_fixed3[i]]
}


### THIS IS THE WAY TO DYNAMICALLY NAME COLUMNS!
# png.neand[, eval(parse(text=sprintf('p1 := hap_%d_window_pval_table', 1)))]
# png.neand[1:10, eval(parse(text=sprintf('list(p1 = hap_%d_window_pval_table)', 2)))]
# png.neand[1:10][eval(parse(text=sprintf('hap_%d_window_match_N_table >= 100',1)))]
# png.neand[1:10][hap_1_window_match_N_table >= 100]

cat('Split by haplotypes 1 and 2...\n')

extract_haps = function(dt) {
  dt.haps = data.table()
  for (hap in c(1,2)) {
    select_txt = sprintf('hap_%d_window_match_N_table >= 100 & hap_%d_window_match_mh_table >= 10 & hap_%d_n_s_star_snps > 2',
                         hap, hap, hap)
    
    col.names = names(dt)
    col.names.haps = col.names[grepl(sprintf('hap_%d', hap), col.names)]
    col.names.haps.generic = sub(sprintf('hap_%d_', hap), 'hap_', col.names.haps)
    col.names.nohaps = col.names[!grepl('hap_', col.names)]
    
    columns_txt = paste( paste(col.names.nohaps, collapse=', '), 
                         sprintf('hap = %d', hap),
                         paste(sprintf('%s = %s', col.names.haps.generic, col.names.haps), collapse=', '), 
                         sep=', ')
    columns_txt = sprintf('list(%s)', columns_txt)
    columns_txt
    
    dt.haps = rbind(dt.haps, 
                    dt[eval(parse(text=select_txt)), eval(parse(text=columns_txt))])
  }
  setkey(dt.haps, chrom, winstart, winend, ind_id, hap)
  dt.haps
}

setnames(png.neand, 'n_s_star_snps_hap1', 'hap_1_n_s_star_snps')
setnames(png.neand, 'n_s_star_snps_hap2', 'hap_2_n_s_star_snps')
png.neand.haps = extract_haps(png.neand)

write.table(png.neand.haps, gzfile(sprintf('%s.processed_haps.gz', neand_file)), sep='\t', quote=F, col.names=T, row.names=F)
