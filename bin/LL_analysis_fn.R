
f <- function(xy, n, x, y, ...) {
  #
  # Estimate the total where the density does not exceed that at (x,y).
  # from: http://stats.stackexchange.com/questions/63447/integrating-kernel-density-estimator-in-2d
  # I am just interested in z (and the image, for viewing results)
  #
  # `xy` is a 2 by ... array of points.
  # `n`  specifies the numbers of rows and columns to use.
  # `x` and `y` are coordinates of "probe" points.
  # `...` is passed on to `kde2d`.
  #
  # Returns a list:
  #   image:    a raster of the kernel density
  #   integral: the estimates at the probe points.
  #   density:  the estimated densities at the probe points.
  #
  xy.kde <- kde2d(xy[1,], xy[2,], n=n, ...)
  xy.im <- im(t(xy.kde$z), xcol=xy.kde$x, yrow=xy.kde$y) # Allows interpolation $
  z <- interp.im(xy.im, x, y)                            # Densities at the probe points
  #c.0 <- sum(xy.kde$z)                                   # Normalization factor $
  #i <- sapply(z, function(a) sum(xy.kde$z[xy.kde$z < a])) / c.0
  return(list(image=xy.im, density=z))
}


getLLcalls <- function(outputdir, keys.arc2, keys.arc3, dt.model, ss_thresh, dt.test, tag, lims, plot_title,
                       h.null = 2, h.alt = 4,
                       xpts.1 = seq(.4,.6,.005), xpts.2 = seq(0,.15,.005)) {
  
  setkey(dt.model, ancNe, mhNe, arcNe, arcDiv)
  setkey(keys.arc2, ancNe, mhNe, arcNe, arcDiv)
  setkey(keys.arc3, ancNe, mhNe, arcNe, arcDiv)
  
  ## get the null model set
  dt.null = with(rbind(dt.model[keys.arc2], dt.model[keys.arc3])[ms_intr_sites.arc2 == 0 & ms_intr_sites.arc3 == 0], 
                 data.table(logit(hap_pval.arc2, min = .00000001, max = .9999999),
                            logit(hap_pval.arc3, min = .00000001, max = .9999999),
                            hap_pval.arc2, hap_pval.arc3))
  dt.null = dt.null[complete.cases(dt.null)]
  
  ## get the alt model set for arc2
  dt.alt.arc2 = with(dt.model[keys.arc2][ms_intr_ss_sites.arc2 > 10], 
                     data.table(logit(hap_pval.arc2, min = .00000001, max = .9999999),
                                logit(hap_pval.arc3, min = .00000001, max = .9999999),
                                hap_pval.arc2, hap_pval.arc3))
  dt.alt.arc2 = dt.alt.arc2[complete.cases(dt.alt.arc2)]
  
  ## get the alt model set for arc3
  ## since we always have introgression from arc2 in the simulations, we just set arc2 and arc3 swapped
  dt.alt.arc3 = with(dt.model[keys.arc3][ms_intr_ss_sites.arc2 > 10], 
                     data.table(logit(hap_pval.arc3, min = .00000001, max = .9999999),
                                logit(hap_pval.arc2, min = .00000001, max = .9999999),
                                hap_pval.arc2 = hap_pval.arc3, hap_pval.arc3 = hap_pval.arc2))
  dt.alt.arc3 = dt.alt.arc3[complete.cases(dt.alt.arc3)]
  
  
  ## get the test set at the particular threshold
  dt.test = dt.test[ss_sig.neand >= ss_thresh,
                    list(logit(hap_window_pval_table.neand, min=.00000001, max=.9999999), 
                         logit(hap_window_pval_table.den, min=.00000001, max=.9999999),
                         hap_window_pval_table.neand,
                         hap_window_pval_table.den,
                         ss_sig.neand, ss_sig.den)]
  dt.test = dt.test[complete.cases(dt.test)]
  
  n=200
  #h.null=2
  #h.alt=10
  
  ## get the kernel fits for each model
  f.null = f(t(as.matrix(dt.null)), n, dt.test[[1]], dt.test[[2]], h=h.null, lims=lims)
  f.alt.1 = f(matrix(c(dt.alt.arc2[[1]], dt.alt.arc2[[2]]), nrow = 2, byrow = T),
              n, dt.test[[1]], dt.test[[2]], h=h.alt, lims=lims)
  f.alt.2 = f(matrix(c(dt.alt.arc3[[1]], dt.alt.arc3[[2]]), nrow = 2, byrow = T),
              n, dt.test[[1]], dt.test[[2]], h=h.alt, lims=lims)
  
  ## the proportions of alt1 and alt2 to calculate the LL
  # xpts.1 = seq(0,.75,.025)
  # xpts.1 = seq(.4,.6,.005)
  xpts.1[xpts.1 == 0] = .00001
  xpts.1[xpts.1 == 1] = 1-.00001
  # xpts.2 = xpts.1
  # xpts.2 = seq(0,.15,.005)
  xpts.2[xpts.2 == 0] = .00001
  xpts.1.rep = rep(xpts.1,each=length(xpts.2))
  xpts.2.rep = rep(xpts.2,length(xpts.1))
  
  ## the likelihood function!  given priors, and the kernel fits of the test data for null and both alts
  lkfn <- function(prior.alt1, prior.alt2, f.null, f.alt.1, f.alt.2) {
    ret = list()
    if (prior.alt1 + prior.alt2 >= 1) {
      ret$LL = NA
      return(ret)
    }
    
    prior.null = 1-prior.alt1-prior.alt2
    p_z.null = f.null$density * prior.null
    p_z.alt1 = f.alt.1$density * prior.alt1
    p_z.alt2 = f.alt.2$density * prior.alt2
    #dt.test[p_z.null == 0 & p_z.alt1 == 0 & p_z.alt2 == 0]
    #dt.test[p_z.null == 0]
    #is.na(p_z.null == 0)
    ret$nohits = p_z.null == 0 | is.na(p_z.null) | p_z.alt1 == 0 | is.na(p_z.alt1) | p_z.alt2 == 0 | is.na(p_z.alt2)
    p_z.null[p_z.null == 0 | is.na(p_z.null)] = 1/nrow(dt.test)
    p_z.alt1[p_z.alt1 == 0 | is.na(p_z.alt1)] = 1/nrow(dt.test)
    p_z.alt2[p_z.alt2 == 0 | is.na(p_z.alt2)] = 1/nrow(dt.test)
    ret$LL = sum(log(p_z.null + p_z.alt1 + p_z.alt2))
    ret$post.p.null = p_z.null / (p_z.null + p_z.alt1 + p_z.alt2)
    ret$post.p.alt1 = p_z.alt1 / (p_z.null + p_z.alt1 + p_z.alt2)
    ret$post.p.alt2 = p_z.alt2 / (p_z.null + p_z.alt1 + p_z.alt2)
    ret$prior.alt1 = prior.alt1
    ret$prior.alt2 = prior.alt2
    ret
  }
  
  ## apply the likelihood function across all of the proportions for alt1 and alt2
  xm = mapply(function(x,y) lkfn(x,y,f.null,f.alt.1,f.alt.2)[['LL']], xpts.1.rep, xpts.2)

  write.table(data.table(LL = xm, prior1 = xpts.1.rep, prior2 = xpts.2), paste0(outputdir, '/LL.callset', tag, '.fullLL.txt'), sep="\t", col.names = T, row.names=F, quote=F)

  save(xm, file=sprintf('%s/xm.%s.Rdata', outputdir, tag))
  xm = matrix(xm, nrow = length(xpts.1))
  
  ## now save some critical info for returning [includes rerunning on the optimal proportions]
  max.ix = which(xm == max(xm, na.rm = T))
  max.ret = lkfn(xpts.1.rep[max.ix], xpts.2.rep[max.ix], f.null, f.alt.1, f.alt.2)
  
  ## make figures if required
  
  pdf(sprintf('%s/test_LL_classifier.null_alt.%s.image.pdf', outputdir, tag))
  image(f.null$image)
  abline(v=-5, h=-5)
  #title(plot_title)
  image(f.alt.1$image)
  image(f.alt.2$image)
  #image(f.alt.1$image + f.alt.2$image)
  abline(v=-5, h=-5)
  #image(f.null$image + f.alt.1$image + f.alt.2$image)
  abline(v=-5, h=-5)
  dev.off()
  
  xm.trim = xm
  xm.trim[xm < quantile(xm, .9, na.rm=T)] = NA
  
  pdf(sprintf('%s/test_LL_classifier.png.%s.pct_neand_and_den.image.ss_%f.pdf', 
              outputdir, tag, ss_thresh))
  image(1:length(xpts.2), 1:length(xpts.1), matrix(xm.trim, nrow = length(xpts.2)), axes=F)
  axis(1, at = 1:ncol(xm.trim), labels=xpts.2)
  axis(2, at = 1:nrow(xm.trim), labels=xpts.1)
  title(sprintf('Max LL: %.0f', 
                max(xm, na.rm = T)))
  #   title(sprintf('mhNe=%d, ancNe=%d \n arcNe=%d, arcDiv=%d, arcNe=%d, arcDiv=%d\nMax LL: %.0f', 
  #                 dt.keys.outer[key.idx.outer, mhNe], dt.keys.outer[key.idx.outer, ancNe], 
  #                 dt.keys.inner[key.idx.inner.arc2, arcNe], dt.keys.inner[key.idx.inner.arc2, arcDiv], 
  #                 dt.keys.inner[key.idx.inner.arc3, arcNe], dt.keys.inner[key.idx.inner.arc3, arcDiv], 
  #                 max(xm, na.rm = T)))
  image(1:length(xpts.2), 1:length(xpts.1), matrix(xm, nrow = length(xpts.2)), axes=F)
  axis(1, at = 1:ncol(xm), labels=xpts.2)
  axis(2, at = 1:nrow(xm), labels=xpts.1)
  dev.off()
  
  
  hist(dt.test[ss_sig.neand >= ss_thresh, hap_window_pval_table.neand])
  qobj.neand = qvalue::qvalue(dt.test[ss_sig.neand >= ss_thresh, hap_window_pval_table.neand])
  names(qobj.neand)
  qobj.neand$pi0
  1-max.ret$prior.alt1-max.ret$prior.alt2
  
  hist(dt.test[ss_sig.neand >= ss_thresh, hap_window_pval_table.den])
  hist(dt.test[ss_sig.den >= ss_thresh, hap_window_pval_table.den])
  qobj.den = qvalue::qvalue(dt.test[ss_sig.neand >= ss_thresh, hap_window_pval_table.den])
  names(qobj.den)
  qobj.den$pi0
  1-max.ret$prior.alt1-max.ret$prior.alt2
  
  #
  
  dt.test[, null.prob := max.ret$post.p.null]
  dt.test[, alt1.prob := max.ret$post.p.alt1]
  dt.test[, alt2.prob := max.ret$post.p.alt2]
  for (pt in c(1, .9, .5, .1)) {
    dt.test[null.prob <= pt, null.cat := pt]
    dt.test[alt1.prob <= pt, alt1.cat := pt]
    dt.test[alt2.prob <= pt, alt2.cat := pt]
  }
  
  
  
  null_thresholds = seq(0,.8,.001)
  fdrs = sapply(null_thresholds, function(null_thresh) dt.test[null.prob < null_thresh, sum(null.prob) / .N])
  anti_fdrs = sapply(null_thresholds, function(null_thresh) dt.test[null.prob > null_thresh, sum(alt1.prob + alt2.prob) / .N])

  pct_non_null = sapply(null_thresholds, function(null_thresh) dt.test[, sum(null.prob < null_thresh) / .N])
  pct_non_null.alt1 = sapply(null_thresholds, function(null_thresh) dt.test[, sum(null.prob < null_thresh & alt1.prob > .5) / .N])
  pct_non_null.alt2 = sapply(null_thresholds, function(null_thresh) dt.test[, sum(null.prob < null_thresh & alt2.prob > .5) / .N])
  pct_non_null.ambig = sapply(null_thresholds, function(null_thresh) dt.test[, sum(null.prob < null_thresh & alt1.prob < .5 & alt2.prob < .5) / .N])

  write.table(data.table(null_thresholds, fdrs, anti_fdrs, pct_non_null, pct_non_null.alt1, pct_non_null.alt2, pct_non_null.ambig), 
              paste0(outputdir, '/LL.callset', tag, '.null_pp_thresh_vs_fdr.txt'), sep="\t", col.names = T, row.names=F, quote=F)

  dt.thresh = data.table(null_thresholds, fdrs, pct_non_null, pct_non_null.alt1, pct_non_null.alt2)
  tail(dt.thresh[fdrs <= .05])
  null_thresholds.05 =  tail(dt.thresh[fdrs <= .05, null_thresholds], n = 1)
  max.ret$null_threshold.05 = null_thresholds.05
    
  pdf(sprintf('%s/test_LL_classifier.png.%s.fdr_thresh_figs.ss_%f.pdf', 
              outputdir, tag, ss_thresh))
  plot(null_thresholds, fdrs)
  abline(h=.05)
  abline(h=.1)
  plot(null_thresholds, pct_non_null)
  plot(fdrs, pct_non_null)
  abline(v=.05)
  abline(v=.1)
  points(fdrs, pct_non_null.alt1, pch=2)
  points(fdrs, pct_non_null.alt2, pch=3)
  points(fdrs, pct_non_null.ambig, pch=3)
  dev.off()
  
  qobj.neand.pthresh = max(qobj.neand$pvalues[qobj.neand$qvalues <= .05])
  sum(qobj.neand$qvalues <= .05) / length(qobj.neand$qvalues)
  qobj.pct_neand = sum(qobj.neand$qvalues <= .05 & qobj.den$qvalues > .05) / length(qobj.neand$qvalues)
  LLobj.pct_neand = dt.test[, sum(alt1.prob > .5 & null.prob < null_thresholds.05) / .N]
  
  qobj.den.pthresh = max(qobj.den$pvalues[qobj.den$qvalues <= .05])
  sum(qobj.den$qvalues <= .05) / length(qobj.den$qvalues)
  qobj.pct_den = sum(qobj.den$qvalues <= .05 & qobj.neand$qvalues > .05) / length(qobj.neand$qvalues)
  dt.test[, sum(alt2.prob > .5) / .N]
  dt.test[, sum(alt2.prob > .5 & null.prob < null_thresholds.05) / .N]
  LLobj.pct_den = dt.test[, sum(alt2.prob > .5 & null.prob < null_thresholds.05) / .N]
  
  max.ret$pct_neand = LLobj.pct_neand
  max.ret$pct_den = LLobj.pct_den
  
  qobj.pct_ambig_arc = sum(qobj.den$qvalues <= .05 & qobj.neand$qvalues <= .05) / length(qobj.neand$qvalues)
  LLobj.pct_ambig_arc = dt.test[, sum(alt2.prob < .5 & alt1.prob < .5 & null.prob < null_thresholds.05) / .N]
  
  cat(null_thresholds.05, '\n')
  cat(sprintf('LL N:%.3f D:%.3f A:%.3f\nqval N:%.3f D:%.3f A:%.3f\n',
              LLobj.pct_neand, LLobj.pct_den, LLobj.pct_ambig_arc,
              qobj.pct_neand, qobj.pct_den, qobj.pct_ambig_arc), '\n')
  
  png(sprintf('%s/test_LL_classifier.png.%s.neand_vs_den_cats.ss_%f.png', 
              outputdir, tag, ss_thresh), width=1000, height=1000)
  # ggplot(dt.test, aes(x=V1, y=V2)) + geom_point(alpha=.05) + facet_wrap(~null.cat)
  print(ggplot(mapping = aes(x=V1, y=V2)) + 
          geom_point(data=dt.test[null.prob >= null_thresholds.05], color='red') +
          geom_point(data=dt.test[null.prob < null_thresholds.05 & alt1.prob/alt2.prob > 2], color='blue') +
          geom_point(data=dt.test[null.prob < null_thresholds.05 & alt2.prob/alt1.prob > 2], color='green') +
          geom_vline(xintercept = logit(qobj.neand.pthresh), linetype=2) +
          geom_hline(yintercept = logit(qobj.den.pthresh), linetype=2) +
          geom_vline(xintercept = -5) +
          geom_hline(yintercept = -5) +
          ggtitle(sprintf('%s\nLL N:%.3f D:%.3f A:%.3f\nqval N:%.3f D:%.3f A:%.3f\n',
                          plot_title, 
                          LLobj.pct_neand, LLobj.pct_den, LLobj.pct_ambig_arc,
                          qobj.pct_neand, qobj.pct_den, qobj.pct_ambig_arc)) +
          xlab('Neanderthal logit p-values') +
          ylab('Denisovan logit p-values'))
    dev.off()
  
  
  
       # pdf(sprintf('%s/%s.png.%d.%d.%d.pct_neand_and_den.persp.ss_%f.pdf', outputdir, 
       #             dt.keys[key.idx, mhNe], dt.keys[key.idx, arcNe], dt.keys[key.idx, ancNe], 
       #             thresh))
       # print(persp(matrix(xm, nrow = length(xpts.2)), theta=-50, phi=20))
       # title(sprintf('mhNe=%d, arcNe=%d, ancNe=%d\nMax LL: %.0f', 
       #               dt.keys[key.idx, mhNe], dt.keys[key.idx, arcNe], dt.keys[key.idx, ancNe], 
       #               max(xm, na.rm = T)))
       # dev.off()
       # pdf(sprintf('%s/%s.png.sim_image.%d.%d.%d.pdf', outputdir, 
       #             dt.keys[key.idx, mhNe], dt.keys[key.idx, arcNe], dt.keys[key.idx, ancNe]))
       # image(f.null$image + f.alt.1$image + f.alt.2$image)
       # image(f.alt.1$image + f.alt.2$image)
       # dev.off()

  save.image(sprintf('%s/LL_fn_workspace.%s.Rdata', outputdir, tag))
  
  return(max.ret)
}
