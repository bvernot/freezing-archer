from __future__ import division
import sys
from numpy import std
import itertools
from collections import defaultdict, Counter

local_debug = False

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    from http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def initialize_analysis(opts):
    if opts.ancestral_bsg == None:
        print "D-statistic estimates require an ancestral sequence!"
        sys.exit(-1)
        pass
    pass


def run_window_analysis(chrom, winstart, winend, snps, opts):

    # (((h1, h2), neand), chimp)
    # ABBA
    # BABA
    
    # specifically, looking for sites where:
    # - chimp is ancestral
    # - neand is derived
    # - sites polymorphic between two pops / inds

    abba = defaultdict(int)
    baba = defaultdict(int)

    d_sum = 0
    count = 0

    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    inregion = opts.regions.amount_in_region(full_chrom, winstart, winend, fail=False)
    
    # for ind_t in opts.target_indices:
    #     for ind_r in opts.reference_indices:
    #         ## considering two individuals, one from target, one from ref
    #         # keep all sites where there's a difference between the individual's genotypes (?)
    #         # or actually keep any site where there's some polymorphism?
    #         ind_snps = [(s['genotypes'][ind_t], s['genotypes'][ind_r], s['arc_genotype']) \
    #                         for s in snps \
    #                         if s['arc_is_derived'] and s['genotypes'][ind_t] - s['genotypes'][ind_r] != 0]
    #         # if s['genotypes'][ind_t] - s['genotypes'][ind_r] != 0 or s['genotypes'][ind_t] == 1]
            
    #         abba = sum(g_t > g_r for (g_t,g_r,_) in ind_snps)
    #         baba = sum(g_r > g_t for (g_t,g_r,_) in ind_snps)
            
    #         d = (abba - baba) / (abba + baba) \
    #             if (abba + baba) != 0 else None

    #         # d_num1 = sum(s['genotypes'][ind_t]/2 * (1 - s['genotypes'][ind_r]/2) * s['arc_genotype']/2 - \
    #         #                  s['genotypes'][ind_r]/2 * (1 - s['genotypes'][ind_t]/2) * s['arc_genotype']/2 \
    #         #                  for s in snps)
    #         # d_den1 = sum(s['genotypes'][ind_t]/2 * (1 - s['genotypes'][ind_r]/2) * s['arc_genotype']/2 + \
    #         #                  s['genotypes'][ind_r]/2 * (1 - s['genotypes'][ind_t]/2) * s['arc_genotype']/2 \
    #         #                  for s in snps)

    #         d_num1 = sum(g_t/2 * (1 - g_r/2) * g_a/2 - \
    #                          g_r/2 * (1 - g_t/2) * g_a/2 \
    #                          for (g_t,g_r,g_a) in ind_snps)
    #         d_den1 = sum(g_t/2 * (1 - g_r/2) * g_a/2 + \
    #                          g_r/2 * (1 - g_t/2) * g_a/2 \
    #                          for (g_t,g_r,g_a) in ind_snps)

    #         print 'DSTATa', chrom, winstart, winend, inregion, ind_t, ind_r, abba, baba, d
    #         print 'DSTATb', chrom, winstart, winend, inregion, ind_t, ind_r, d_num1, d_den1, d_num1/d_den1 if d_den1 != 0 else None
            
    #         if d is not None:
    #             d_sum += d
    #             count += 1
    #             pass

    #         pass
    #     pass

    # print 'DAVG', chrom, winstart, winend, inregion, len(snps), count, d_sum/count if count != 0 else 'NA'


    # for i, tinds in enumerate(chunks(opts.target_indices, 20)):
    #     for j, rinds in enumerate(chunks(opts.reference_indices, 20)):
    #         if len(tinds) != 20 or len(rinds) != 20: continue
            
    #         ind_snps = [(sum(s['genotypes'][ind_t] for ind_t in tinds), \
    #                          sum(s['genotypes'][ind_r] for ind_r in rinds), \
    #                          s['arc_genotype']) \
    #                         for s in snps \
    #                         if s['arc_is_derived']]

    #         #abba = sum(g_t > 0 and g_r == 0 for (g_t,g_r,_) in ind_snps)
    #         #baba = sum(g_r > 0 and g_t == 0 for (g_t,g_r,_) in ind_snps)
    #         abba = sum(g_t > g_r for (g_t,g_r,_) in ind_snps)
    #         baba = sum(g_r > g_t for (g_t,g_r,_) in ind_snps)
            
    #         d = (abba - baba) / (abba + baba) \
    #             if (abba + baba) != 0 else None

    #         # d_num1 = sum(s['genotypes'][ind_t]/2 * (1 - s['genotypes'][ind_r]/2) * s['arc_genotype']/2 - \
    #         #                  s['genotypes'][ind_r]/2 * (1 - s['genotypes'][ind_t]/2) * s['arc_genotype']/2 \
    #         #                  for s in snps)
    #         # d_den1 = sum(s['genotypes'][ind_t]/2 * (1 - s['genotypes'][ind_r]/2) * s['arc_genotype']/2 + \
    #         #                  s['genotypes'][ind_r]/2 * (1 - s['genotypes'][ind_t]/2) * s['arc_genotype']/2 \
    #         #                  for s in snps)

    #         d_num1 = sum(g_t/40 * (1 - g_r/40) * g_a/2 - \
    #                          g_r/40 * (1 - g_t/40) * g_a/2 \
    #                          for (g_t,g_r,g_a) in ind_snps)
    #         d_den1 = sum(g_t/40 * (1 - g_r/40) * g_a/2 + \
    #                          g_r/40 * (1 - g_t/40) * g_a/2 \
    #                          for (g_t,g_r,g_a) in ind_snps)

    #         print 'D20STATa', chrom, winstart, winend, inregion, i, j, abba, baba, d
    #         print 'D20STATb', chrom, winstart, winend, inregion, i, j, d_num1, d_den1, d_num1/d_den1 if d_den1 != 0 else None
            
    #         pass
    #     pass
        
    

    tinds = opts.target_indices
    rinds = opts.reference_indices

            
    ind_snps = [(sum(s['genotypes'][ind_t] for ind_t in tinds), \
                     sum(s['genotypes'][ind_r] for ind_r in rinds), \
                     s['arc_genotype']) \
                    for s in snps \
                    if s['arc_is_derived']]
    
    d_num1 = sum(g_t/len(tinds)/2 * (1 - g_r/len(rinds)/2) * g_a/2 - \
                     g_r/len(rinds)/2 * (1 - g_t/len(tinds)/2) * g_a/2 \
                     for (g_t,g_r,g_a) in ind_snps)
    d_num1_hom = sum(g_a/2 * (1 - g_r/len(rinds)/2) * g_a/2 - \
                         g_r/len(rinds)/2 * (1 - g_a/2) * g_a/2 \
                         for (g_t,g_r,g_a) in ind_snps)
    d_num1_hom_corrected = sum(g_a/2 if g_a/2 > g_t/len(tinds)/2 else g_t/len(tinds)/2 \
                                   * (1 - g_r/len(rinds)/2) \
                                   * g_a/2 if g_a/2 > g_t/len(tinds)/2 else g_t/len(tinds)/2  - \
                                   g_r/len(rinds)/2 \
                                   * (1 - g_a/2 if g_a/2 > g_t/len(tinds)/2 else g_t/len(tinds)/2) \
                                   * g_a/2 if g_a/2 > g_t/len(tinds)/2 else g_t/len(tinds)/2 \
                                   for (g_t,g_r,g_a) in ind_snps)
    d_den1 = sum(g_t/len(tinds)/2 * (1 - g_r/len(rinds)/2) * g_a/2 + \
                     g_r/len(rinds)/2 * (1 - g_t/len(tinds)/2) * g_a/2 \
                     for (g_t,g_r,g_a) in ind_snps)

    if opts.first_line:
        print 'stat_name', 'chrom', 'start', 'end', 'bases', 'sites', 'd_numer', 'd_denom', 'dstat', 'subpop', 'refpop', 'ninds_target', 'ninds_ref', 'd_numer_hom', 'd_numer_hom_corrected', \
            ' '.join(opts.tag_ids)
        opts.first_line = False
        pass
        
    # print ind_snps
    print 'DSTATc', chrom, winstart, winend, inregion, len(ind_snps), d_num1, d_den1, d_num1/d_den1 if d_den1 != 0 else None, \
        '.'.join(opts.target_populations), '.'.join(opts.reference_populations), len(opts.target_indices), len(opts.reference_indices), \
        d_num1_hom, d_num1_hom_corrected, ' '.join(opts.tags)

            
    
    return


def finish_analysis(opts):
    return
