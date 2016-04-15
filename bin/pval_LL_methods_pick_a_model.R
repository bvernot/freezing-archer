library('ggplot2')
library('data.table')
#source('~/Dropbox/archaic exomes/multiplot.R')
library('gtools')
library(MASS)
library(spatstat) # im class

options(bitmapType='cairo')

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)


if (length(args) != 4) {
   cat("Requires:\n - population-id (e.g., PNG)\n - ...neand.ALLCHRS.gz.processed_haps.gz\n - ...den.ALLCHRS.gz.processed_haps.gz\n - outputdir\n")
   quit()
}


#mypop = 'SAS'
mypop = args[1]
neand_pvals_file = args[2]
den_pvals_file   = args[3]
outputdir   = args[4]

cat(mypop)


####################
#### REAL PVALS
####################

# dt.pvals.neand = fread(paste0('gunzip -c /net/akey/vol2/bvernot/archaic_1kg_p3/experiments/pvalue_tables/query_trial_output6.new_filters.tsvcat.ALLCHRS',
#                                   '/query_trial_output6.new_filters.tsvcat.', mypop, '.neand.ALLCHRS.gz.processed_haps.gz'))
# dt.pvals.den   = fread(paste0('gunzip -c /net/akey/vol2/bvernot/archaic_1kg_p3/experiments/pvalue_tables/query_trial_output6.new_filters.tsvcat.ALLCHRS',
#                                   '/query_trial_output6.new_filters.tsvcat.', mypop, '.den.ALLCHRS.gz.processed_haps.gz'))

dt.pvals.neand = fread(paste0('gunzip -c ', neand_pvals_file))
dt.pvals.den   = fread(paste0('gunzip -c ', den_pvals_file))

# dt.pvals.neand = fread(paste0('gunzip -c /net/akey/vol2/bvernot/archaic_1kg_p3/experiments/pvalue_tables/query_trial_output6.new_filters.tsvcat.ALLCHRS',
#                                   '/query_trial_output6.new_filters.tsvcat.', mypop, '.neand.ALLCHRS.gz.processed_haps.gz.ss_95.gz'))
# dt.pvals.den   = fread(paste0('gunzip -c /net/akey/vol2/bvernot/archaic_1kg_p3/experiments/pvalue_tables/query_trial_output6.new_filters.tsvcat.ALLCHRS',
#                                   '/query_trial_output6.new_filters.tsvcat.', mypop, '.den.ALLCHRS.gz.processed_haps.gz.ss_95.gz'))

setkey(dt.pvals.neand, chrom, winstart, winend, ind_id, hap)
setkey(dt.pvals.den, chrom, winstart, winend, ind_id, hap)
dt.pvals.merged = merge(dt.pvals.den, dt.pvals.neand, suffixes = c('.den', '.neand'))


####################
#### SIMULATED PVALS
####################
# dt.intr23.arc2 = fread(paste0('/net/akey/vol2/bvernot/archaic_1kg_p3/bin/tsvcat /net/akey/vol2/bvernot/archaic_1kg_p3/experiments/pvalue_tables_from_simulations/',
#                               'pval_output.different_ast.with_mig_tables.summary.iter1-8_noshuf/*.arc_2.intr23.txt'))[, 
#                               .SD, keyby=list(chrom, pop, ast1, ast2, arcNe, ancNe, mhNe, arcDiv, iter, seg, hap)]
# dt.intr23.arc3 = fread(paste0('/net/akey/vol2/bvernot/archaic_1kg_p3/bin/tsvcat /net/akey/vol2/bvernot/archaic_1kg_p3/experiments/pvalue_tables_from_simulations/',
#                               'pval_output.different_ast.with_mig_tables.summary.iter1-8_noshuf/*.arc_3.intr23.txt'))[, 
#                               .SD, keyby=list(chrom, pop, ast1, ast2, arcNe, ancNe, mhNe, arcDiv, iter, seg, hap)]

dt.intr23.arc2 = fread(paste0('bin/tsvcat program_data/',
                              'pval_output.different_ast.with_mig_tables.summary.iter1-8_noshuf/*.arc_2.intr23.txt'))[, 
                              .SD, keyby=list(chrom, pop, ast1, ast2, arcNe, ancNe, mhNe, arcDiv, iter, seg, hap)]
dt.intr23.arc3 = fread(paste0('bin/tsvcat program_data/',
                              'pval_output.different_ast.with_mig_tables.summary.iter1-8_noshuf/*.arc_3.intr23.txt'))[, 
                              .SD, keyby=list(chrom, pop, ast1, ast2, arcNe, ancNe, mhNe, arcDiv, iter, seg, hap)]

dt.intr23.merge = merge(dt.intr23.arc2, dt.intr23.arc3, suffixes = c('.arc2', '.arc3'))




###################
#### MLE functions

source('bin/LL_analysis_fn.R')

#### GET THE VARIOUS COMBOS THAT ARE POSSIBLE (I.E., ARC2 DIVERGED 200KYA, ARC3 DIVERGED 300KYA)

setkey(dt.intr23.merge, ancNe, mhNe)
dt.keys.outer = dt.intr23.merge[, .N, keyby=list(ancNe, mhNe)]

dt.allkeys = data.table()
for (key.idx.outer in seq(1,nrow(dt.keys.outer))) {
  
  dt.outer = dt.intr23.merge[dt.keys.outer[key.idx.outer]]
  setkey(dt.outer, arcNe, arcDiv)
  dt.keys.inner = dt.outer[, .N, keyby=list(arcNe, arcDiv)]
  
  for (key.idx.inner.arc2 in seq(1,nrow(dt.keys.inner))) {
    for (key.idx.inner.arc3 in seq(1,nrow(dt.keys.inner))) {
      
      idx.tag = sprintf('%d.%d_%d.%d.%d.%d',
                        dt.keys.outer[key.idx.outer, mhNe], dt.keys.outer[key.idx.outer, ancNe], 
                        dt.keys.inner[key.idx.inner.arc2, arcNe], dt.keys.inner[key.idx.inner.arc2, arcDiv], 
                        dt.keys.inner[key.idx.inner.arc3, arcNe], dt.keys.inner[key.idx.inner.arc3, arcDiv])
      
      dt.allkeys = rbind(dt.allkeys,
                         cbind(dt.keys.outer[key.idx.outer, list(ancNe, mhNe)], 
                               dt.keys.inner[key.idx.inner.arc2, list(arcNe.1 = arcNe, arcDiv.1 = arcDiv, N1 = N)],
                               dt.keys.inner[key.idx.inner.arc3, list(arcNe.2 = arcNe, arcDiv.2 = arcDiv, N2 = N)],
                               data.table(tag = idx.tag)))
    }
  }
}
dt.allkeys.arc2 = dt.allkeys[, list(ancNe, mhNe, arcNe = arcNe.1, arcDiv = arcDiv.1)]
dt.allkeys.arc3 = dt.allkeys[, list(ancNe, mhNe, arcNe = arcNe.2, arcDiv = arcDiv.2)]
dt.allkeys



################
#### SET UP PARAMS FOR LL FUNCTION

lims <- c(min(dt.intr23.merge[, logit(hap_pval.arc2, min = .00000001, max = .9999999)],
              dt.intr23.merge[, logit(hap_pval.arc3, min = .00000001, max = .9999999)], na.rm=T)-5, 
          max(dt.intr23.merge[, logit(hap_pval.arc2, min = .00000001, max = .9999999)],
              dt.intr23.merge[, logit(hap_pval.arc3, min = .00000001, max = .9999999)], na.rm=T)+5)
lims = c(lims, lims)


dt.LL = data.table()
thresh=0.9900
h.null=1.5
h.alt=6

# keyline = dt.allkeys[, which(arcDiv.1 == 200000 & arcDiv.2 == 350000)]
# keyline = dt.allkeys[, which(arcDiv.1 == 200000 & arcDiv.2 == 300000)]

# keylines = c(dt.allkeys[, which(arcDiv.1 == 200000 & arcDiv.2 == 350000)],
#              dt.allkeys[, which(arcDiv.1 == 200000 & arcDiv.2 == 300000)],
#              dt.allkeys[, which(arcDiv.1 == 150000 & arcDiv.2 == 350000)],
#              dt.allkeys[, which(arcDiv.1 == 150000 & arcDiv.2 == 300000)])
# keylines = seq(1, nrow(dt.allkeys))

######
### for PNG sequences, the best model fit was with N divergence = 150kya, and D div = 350kya
keyline = dt.allkeys[, which(arcDiv.1 == 150000 & arcDiv.2 == 350000)]

#for (keyline in keylines) {
for (thresh in c(0.9900, 0.9800, 0.9850, 0.9900, 0.9950)) {
#for (thresh in c(0.9500, 0.9600, 0.9700, 0.9800, 0.9850, 0.9900, 0.9950)) {

      keys.arc2 = dt.allkeys.arc2[keyline]
      keys.arc3 = dt.allkeys.arc3[keyline]

      ######
      ### calculate LL for different proportions N and D in the S* callset
      ### these lines search a very coarse grid (increments of 5%)
      xpts.1 = seq(0,.65,.005)
      xpts.2 = seq(0,.25,.005)
      # xpts.1 = seq(0,.65,.2)
      # xpts.2 = seq(0,.65,.2)
      ## you can search a much finer grid, but it will take a long time!
      # xpts.1 = seq(0,.65,.005)
      # xpts.2 = seq(0,.65,.005)

      tag = paste0(mypop, '_LL_', thresh, '_iter.90_largetables_fulldata_restricted_range_', 'null', h.null, '_', 'alt', h.alt, '_', dt.allkeys[keyline]$tag)

            mr.sam = getLLcalls(outputdir, keys.arc2, keys.arc3, dt.intr23.merge, thresh, 
                                dt.pvals.merged,
                                tag, lims,
                                paste0(mypop, ' Archaic Calls; S* thresh ', thresh, '; mig tables\n', dt.allkeys[keyline]$tag),
                                h.null = h.null, h.alt = h.alt,
                                xpts.1 = xpts.1,
                                xpts.2 = xpts.2)
            dt.LL = rbind(dt.LL, data.table(LL = mr.sam$LL, prior.alt1 = mr.sam$prior.alt1, 
                                            prior.alt2 = mr.sam$prior.alt2, dt.allkeys[keyline], thresh,
                                            h.null, h.alt, pop = mypop))
          
            write.table(dt.LL, paste0(outputdir, '/LL.callset', tag, '.txt'), sep="\t")

      save(mr.sam, file=paste0(outputdir, '/LL.callset', tag, '.Rdata'))
      # load(paste0('~/Dropbox/LL.callset', mypop, '.mr_', thresh, '.Rdata'))
      # load(paste0('~/Dropbox/LL.callset.mr_', thresh, '.Rdata'))
      dt.test = dt.pvals.merged[ss_sig.neand >= thresh,
                                       list(logit(hap_window_pval_table.neand, min=.00000001, max=.9999999), 
                                            logit(hap_window_pval_table.den, min=.00000001, max=.9999999),
                                            hap_window_pval_table.neand,
                                            hap_window_pval_table.den,
                                            ind_id, hap,
                                            ss_sig.neand, ss_sig.den, chrom, winstart, winend,
                                            hap_s_start.neand, hap_s_end.neand,
                                            hap_s_start.den, hap_s_end.den)]
      dt.test[, in_test_set := complete.cases(dt.test)]
      dt.test[, .N]
      dt.test[in_test_set == F]
      dt.test[in_test_set == T, post.p.null := mr.sam$post.p.null]
      dt.test[in_test_set == T, post.p.alt1 := mr.sam$post.p.alt1]
      dt.test[in_test_set == T, post.p.alt2 := mr.sam$post.p.alt2]
      # pdf(file=paste0(outputdir, '/LL.callset', mypop, '.mr_chr1_', thresh, '.pdf'))
      # print(ggplot(dt.test[chrom == 1 & post.p.null <= mr.sam$null_threshold.05], aes(x=winstart)) + geom_histogram(binwidth=1e6) +
      #   geom_vline(xintercept=c(110e6, 125e6)))
      # dev.off()

      write.table(dt.test[post.p.null <= mr.sam$null_threshold.05,
                             list(chrom, winstart, winend, ind_id, hap,
                                  sprintf('%d.%s.%d', chrom, ind_id, hap),
                                  hap_s_start.den, hap_s_end.den,
                                  hap_s_start.neand, hap_s_end.neand,
                                  ss_sig.neand, ss_sig.den,
                                  hap_window_pval_table.neand,
                                  hap_window_pval_table.den,
                                  post.p.null, post.p.alt1, post.p.alt2)],
                  file=paste0(outputdir, '/LL.callset', tag, '.bed'),
                  sep = '\t', row.names = F, quote=F)

      write.table(dt.test[,
                             list(chrom, winstart, winend, ind_id, hap,
                                  sprintf('%d.%s.%d', chrom, ind_id, hap),
                                  hap_s_start.den, hap_s_end.den,
                                  hap_s_start.neand, hap_s_end.neand,
                                  ss_sig.neand, ss_sig.den,
                                  hap_window_pval_table.neand,
                                  hap_window_pval_table.den,
                                  post.p.null, post.p.alt1, post.p.alt2)],
                  file=paste0(outputdir, '/LL.callset', tag, '.FULL.bed'),
                  sep = '\t', row.names = F, quote=F)


}
