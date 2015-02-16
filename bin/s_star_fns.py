from __future__ import division
import sys
from numpy import std
import itertools
from collections import defaultdict, Counter

s_star_debug = False

def calc_geno_dist(gt1, gt2):
    gd = [abs(gt1[i] - gt2[i]) for i in range(len(gt1))]
    return sum(gd)

def calc_s(gt1, gt2, p1, p2):

    if abs(p2-p1) < 10:
        return -sys.maxint

    gd = calc_geno_dist(gt1, gt2)
    if gd == 0:
        return 5000 + abs(p2-p1)
    if gd < 5:
        return -10000
    
    return -sys.maxint

def calc_s_star(genotypes, positions, nsnps):
    
    s_star_scores = [0] * nsnps
    s_star_snps = [[]] * nsnps

    ## loop through all snps, calculating s_star given that:
    #### you have a set of snps that ends in k
    #### you are now trying to figure out the set of snps before k that maximizes the score with k
    ## try all possible j (i.e., sets that end in j), save the best one
    ## consider both the previous best set with j, and that you start over with (j,k)

    for k in range(nsnps):

        max_score = -sys.maxint
        max_snps = []
        for j in range(k):
            # just the set (j,k)
            score1 = calc_s(genotypes[j], genotypes[k], positions[j], positions[k])
            # the previous set that ends in j, plus k
            score2 = s_star_scores[j] + score1

            if s_star_debug: print 'current snps (j=%d,k=%d):' % (j, k), 'max_score=%d, score for %s = %d, score for %s = %d' % (max_score, (j,k), score1, s_star_snps[j] + [k], score2)
            if s_star_debug: print '  best score for j:', s_star_scores[j]
            if s_star_debug: print '  best snps  for j:', s_star_snps[j]

            if max_score < score2:
                max_score = score2
                max_snps = s_star_snps[j] + [k]
                if s_star_debug: print '   settting max_score, snps, keeping j set       :', max_score, max_snps
                pass

            if max_score < score1:
                max_score = score1
                max_snps = [j, k]
                if s_star_debug: print '   settting max_score, snps, starting over with j:', max_score, max_snps
                pass
            pass

        s_star_scores[k] = max_score
        s_star_snps[k] = max_snps

        if s_star_debug: print
        if s_star_debug: print 'current best scores (k=%d)' % k, s_star_scores
        if s_star_debug: print 'current best snps   (k=%d)' % k, s_star_snps
        if s_star_debug: print

        pass
    
    s_star = max(s_star_scores)
    chosen_end_snp = s_star_scores.index(s_star)
    chosen_snps = s_star_snps[chosen_end_snp]
    
    return( s_star, chosen_snps )


def calc_match_pval_from_genome_tables(chrom, winstart, winend, ind_snps, opts):

    match_pval_debug = False
    
    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    my_mapped_bases = opts.regions.amount_in_region(full_chrom, winstart, winend) if opts.regions != None else opts.window_length
    my_mapped_bases_bin = my_mapped_bases // 1000 * 1000

    # get neand hom snps for this region - i.e., we're only looking for sites that have only non-ref
    neand_pos = [p for p in range(winstart, winend) \
                     if opts.neand_vcf.has_var_at_site(chrom, p) and \
                     not opts.neand_vcf.has_ref(chrom, p) and \
                     (opts.regions == None or opts.regions.in_region_one_based(full_chrom, p))]

    # get the number of sites in the region (for this ind, and for neand)
    ind_pos = set([s['pos'] for s in ind_snps])
    ind_pos2 = ind_pos.union(set(neand_pos))
    my_n_sites = len(ind_pos2)

    # count the number of matches to neanderthal
    n_match = sum([s['arc_match'] for s in ind_snps])
    test_ratio = n_match / my_n_sites
    
    # get the average freq
    my_sfs = int(round(sum([s['sfs_target'] for s in ind_snps]) / len(ind_snps)))
    
    # get the std freq
    my_std_dev = int(round(std([s['sfs_target'] for s in ind_snps])))

    # get the distribution of "null" matches

    # first get hits, progressively building up to more and more lax restrictions?
    idx_factor = 2

    # cycle through files
    # for b in set(x for x in (my_mapped_bases_bin - 1000 * idx_factor, my_mapped_bases_bin + 1000 * idx_factor) if x in tbl.keys()):
    if match_pval_debug: print my_mapped_bases_bin, idx_factor
    my_counter = Counter()
    my_null = []
    for idx_factor in range(5):
        for tbl_n, tbl in enumerate(opts.match_pval_table):
            if match_pval_debug: print 'table', tbl_n
            for b in [x for x in range(my_mapped_bases_bin - 2000*idx_factor, my_mapped_bases_bin + 2000*idx_factor+1, 1000) if x in tbl.keys()]:
                if match_pval_debug: print b
                for p in [x for x in range(my_n_sites - 3*idx_factor, my_n_sites + 3*idx_factor+1) if x in tbl[b].keys()]:
                    if match_pval_debug: print b,p
                    for sfs in [x for x in range(my_sfs - 1*idx_factor, my_sfs + 1*idx_factor+1) if x in tbl[b][p].keys()]:
                        if match_pval_debug: print b,p,sfs
                        for sd in [x for x in range(my_std_dev - 5*idx_factor, my_std_dev + 5*idx_factor+1) if x in tbl[b][p][sfs].keys()]:
                            if match_pval_debug: print b,p,sfs,sd,tbl[b][p][sfs][sd]
                            for nm,c in tbl[b][p][sfs][sd].items():
                                my_null += [nm/p] * c
                                if match_pval_debug: print b,p,sfs,sd,nm,c,my_null
                                pass
                            my_counter += tbl[b][p][sfs][sd]
                            pass
                        pass
                    pass
                pass
            if match_pval_debug: print "current counter:", my_counter
            if match_pval_debug: print "current null:", len(my_null), my_null
            pass
        if match_pval_debug: print "finished all files with idx_factor=%d.." % idx_factor
        if len(my_null) > 100: 
            break
        pass

    if match_pval_debug: print "FINISHED GENERATING NULL"
    for n in sorted(my_null):
        if match_pval_debug: print n
        pass

    null_hits = sum([n >= test_ratio for n in my_null])
    total_hits = len(my_null)
    pass

    if total_hits == 0:
        pval = 'NA'
    else:
        pval = null_hits / total_hits
        pass

    # 34396   42      92      63      n_match?  4254    4254    1.0     1
    return (my_mapped_bases, my_n_sites, my_sfs, my_std_dev, n_match, len(ind_snps), total_hits, null_hits, pval, idx_factor)


local_debug = False

def initialize_analysis(opts):
    pass

def calc_local_match_pval(chrom, winstart, winend, snps, opts, ind_indices = None, subset_start = None, subset_end = None):

    if ind_indices == None: ind_indices = opts.target_indices
    if isinstance(ind_indices, int): ind_indices = [ind_indices]

    if subset_start == None: subset_start = winstart
    # subset_end is *inclusive*, unlike winend
    if subset_end == None: subset_end = winend-1

    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    # if we use mapped bases, have to account for subset..
    #my_mapped_bases = opts.regions.amount_in_region(full_chrom, winstart, winend) if opts.regions != None else opts.window_length
    #my_mapped_bases_bin = my_mapped_bases // 1000 * 1000

    # get neand hom snps for this region - i.e., we're only looking for sites that have only non-ref
    neand_pos = [p for p in range(subset_start, subset_end+1) \
                     if opts.neand_vcf.has_var_at_site(chrom, p) and \
                     not opts.neand_vcf.has_ref(chrom, p) and \
                     (opts.regions == None or opts.regions.in_region_one_based(full_chrom, p))]

    # print neand_pos
    # for p in [2000, 2001, 2309, 7879, 16249]:
    #     print '%6d ref=%r var=%r geno=%s' % (p, opts.neand_vcf.has_ref(chrom, p), opts.neand_vcf.has_var_at_site(chrom, p))
    #     pass

    # save two spots for each individual (one for each haplotype)
    # the code below assumes that the ref inds are first, then target (where we actually calculate the pval, not the match pct loop immediately below)
    # but it would be good to factor that out - to just use a fn that returns the indices for target set, for ref set, etc
    match_pct = [None,None] * opts.num_samples

    ## just look at sites in the window of interest
    subset_snps = [s for s in snps if s['pos'] >= subset_start and s['pos'] <= subset_end]
    # print 'all snps %6d %6d' % (winstart, winend), [s['pos'] for s in snps]
    # print 'sub snps %6d %6d' % (subset_start, subset_end), [s['pos'] for s in subset_snps]
    # # print 'sub difs %6d %6d' % (subset_start, subset_end), [subset_end-s['pos'] for s in subset_snps]

    ## loop through all individuals (target and ref) and calculate match pct
    for ind in ind_indices + opts.reference_indices:

        # get haplotype snps (two haplotypes per person, natch)
        ind_hap1 = [s for s in subset_snps if s['haplotypes_1'][ind] > 0]
        ind_hap2 = [s for s in subset_snps if s['haplotypes_2'][ind] > 0]

        # get the number of sites to consider for this haplotype (i.e., on this haplotype and/or in neand)
        ind_pos1 = set([s['pos'] for s in ind_hap1])
        ind_pos1_n = ind_pos1.union(set(neand_pos))
        my_n_sites1 = len(ind_pos1_n)

        ind_pos2 = set([s['pos'] for s in ind_hap2])
        ind_pos2_n = ind_pos2.union(set(neand_pos))
        my_n_sites2 = len(ind_pos2_n)


        # count the number of matches to neanderthal
        n_match1 = sum([s['arc_match'] for s in ind_hap1])
        n_match2 = sum([s['arc_match'] for s in ind_hap2])
        
        # get the pct of sites that match neand
        test_ratio1 = None if my_n_sites1 == 0 else n_match1 / my_n_sites1
        test_ratio2 = None if my_n_sites2 == 0 else n_match2 / my_n_sites2

        # save those values
        match_pct[ind*2] = test_ratio1
        match_pct[ind*2+1] = test_ratio2

        if opts.debug: print 'debug_match_pct_fn IND=%d HAP=1 PCT=%f MATCH_NEAND=%s NEAND_POS=%s DENOM=%s' % (ind, test_ratio1,
                                                                                               str([s['pos'] for s in ind_hap1 if s['arc_match']]),
                                                                                               str(neand_pos),
                                                                                               str(my_n_sites1))
        if opts.debug: print 'debug_match_pct_fn IND=%d HAP=2 PCT=%f MATCH_NEAND=%s NEAND_POS=%s DENOM=%s' % (ind, test_ratio2,
                                                                                               str([s['pos'] for s in ind_hap2 if s['arc_match']]),
                                                                                               str(neand_pos),
                                                                                               str(my_n_sites2))

        if opts.debug: print "test_ratio1", winstart, winend, ind, opts.get_pop_from_sample_index(ind), n_match1, len(neand_pos), len(ind_pos1), my_n_sites1, test_ratio1
        if opts.debug: print "test_ratio2", winstart, winend, ind, opts.get_pop_from_sample_index(ind), n_match2, len(neand_pos), len(ind_pos2), my_n_sites2, test_ratio2

        pass


    # THIS SHOULD BE CONVERTED TO WORK WITH OPTS.TARGET_INDICES (AND MAYBE I SHOULD MAKE A SEPARATE .TARGET_INDICES_HAPS, SO I DON'T NEED TO DO I*2, I*2+1?)
    ref_nones = sum(p == None for p in match_pct[opts.num_target*2:])

    # make a longer pval list than you need, so you can use standard indices
    pvals = [None,None] * opts.num_samples
    # these aren't in the correct order - i.e., one individual's two haplotypes aren't paired with each other
    # but that shouldn't matter - we're always just looking at them in aggregate
    ref_match_pct = [match_pct[i*2] for i in opts.reference_indices] + [match_pct[i*2+1] for i in opts.reference_indices]
    for ind in ind_indices:
        ps_1 = sum(match_pct[ind*2] <= p for p in ref_match_pct)
        ps_2 = sum(match_pct[ind*2+1] <= p for p in ref_match_pct)
        pvals[ind*2] = (ps_1 if ps_1 > 0 else 1) / opts.num_reference / 2
        pvals[ind*2+1] = (ps_2 if ps_2 > 0 else 1) / opts.num_reference / 2
        if opts.debug: print "pval1", winstart, winend, ind, opts.get_pop_from_sample_index(ind), \
                match_pct[ind*2], pvals[ind*2]
        if opts.debug: print "pval2", winstart, winend, ind, opts.get_pop_from_sample_index(ind), \
                match_pct[ind*2+1], pvals[ind*2+1]
        pass

    return (pvals, match_pct, ref_nones)

def run_window_analysis(chrom, winstart, winend, snps, opts):

    if opts.debug: print 'starting window', chrom, winstart, winend, 'NUM SNPS:', len(snps)
    if len(snps) <= 2: return

    (match_pvals2, match_pct2, ref_nones) = calc_local_match_pval(chrom, winstart, winend, snps, opts)

    # # print snps[0]['pos'], snps[-1]['pos']
    # (match_pvals3, match_pct3, ref_nones3) = calc_local_match_pval(chrom, winstart, winend, snps, opts, 0, subset_start = 2300, subset_end = 48988)
    # print 'debug match_pval -  whole region', match_pvals2
    # print 'debug match_pval - subset region', match_pvals3

    # print 'debug match_pct -  whole region', match_pct2
    # print 'debug match_pct - subset region', match_pct3

    # (match_pvals3, match_pct3, ref_nones3) = calc_local_match_pval(chrom, winstart, winend, snps, opts, opts.target_indices)
    # print 'debug match_pval - subset region', match_pvals3
    # print 'debug match_pct - subset region', match_pct3
    

    ## JUST LOOP OVER TARGET INDS
    for ind in opts.target_indices:
        # print snps
        ind1_snps = [s for s in snps if s['genotypes'][ind] > 0 and s['target'] and not s['reference']]
        ind1_snps_all = [s for s in snps if s['genotypes'][ind] > 0]
        ind1_pos = [s['pos'] for s in ind1_snps]

        # ind1_hap1 and ind1_hap2 are just masks for which snps are on which haplotypes - they're the same length as ind1_snps
        ind1_hap1 = [s['haplotypes_1'][ind] > 0 for s in ind1_snps]
        ind1_hap2 = [s['haplotypes_2'][ind] > 0 for s in ind1_snps]


        if opts.debug or local_debug: print "genotypes", [s['genotypes'][ind] for s in ind1_snps]

        if len(ind1_snps) <= 2:
            if local_debug: print 'skipping ind', ind, opts.get_pop_from_sample_index(ind)
            (s_star, s_star_snps) = (0, [])
        else:
            ## s_star is the score, and s_star_snps are the indices (wrt ind1_snps) that are selected by s_star
            (s_star, s_star_snps) = calc_s_star([[s['genotypes'][ind]] for s in ind1_snps], ind1_pos, len(ind1_snps))
            pass

        # for s in s_star_snps:
        #     print s, ind1_snps[s]['pos'], ind1_snps[s]['genotypes'][ind], ind1_snps[s]['arc_match']
        #     pass

        ## find the best haplotype
        n_haps1 = sum(ind1_hap1[i] for i in s_star_snps)
        n_haps2 = sum(ind1_hap2[i] for i in s_star_snps)
        all_haps = (ind1_hap1[i] + ind1_hap2[i]*2 for i in s_star_snps)

        ## currently depreciated..
        if opts.match_pval_table != None:
            match_pval = calc_match_pval(chrom, winstart, winend, ind1_snps_all, opts)
        else:
            match_pval = (None,)
            pass

        if len(s_star_snps) >= 2:
            s_start = ind1_snps[min(s_star_snps)]['pos']
            s_stop  = ind1_snps[max(s_star_snps)]['pos']
            (match_pvals3, match_pct3, ref_nones3) = calc_local_match_pval(chrom, winstart, winend, snps, opts, ind, subset_start = s_start, subset_end = s_stop)
            # print ind, s_start, s_stop, match_pvals3
            # print ind, match_pct3
        else:
            match_pvals3 = ['NA', 'NA'] * opts.num_samples
            match_pct3   = ['NA', 'NA'] * opts.num_samples
            pass

        if opts.first_line:
            print '\t'.join(('chrom', 'winstart', 'winend', 'n_snps', 'n_ind_snps',
                             'ind_id',
                             'pop',
                             's_star',
                             'num_s_star_snps',
                             's_star_snps',
                             'hap_1_window_pval', 'hap_2_window_pval',
                             'hap_1_window_match_pct', 'hap_2_window_match_pct',
                             'hap_1_window_pval_local', 'hap_2_window_pval_local',
                             'hap_1_window_match_pct_local', 'hap_2_window_match_pct_local',
                             'ref_nocals',
                             'n_s_star_snps_hap1', 'n_s_star_snps_hap2',
                             's_star_haps'))
            opts.first_line = False
            pass

        print '\t'.join(str(s) for s in (chrom, winstart, winend, len(snps), len(ind1_snps),
                                         opts.get_id_from_sample_index(ind),
                                         opts.get_pop_from_sample_index(ind),
                                         s_star, 
                                         len(s_star_snps),
                                         ','.join(str(ind1_pos[i]) for i in s_star_snps) if len(s_star_snps) > 0 else '.',
                                         match_pvals2[ind*2], match_pvals2[ind*2+1],
                                         match_pct2[ind*2], match_pct2[ind*2+1],
                                         match_pvals3[ind*2], match_pvals3[ind*2+1],
                                         match_pct3[ind*2], match_pct3[ind*2+1],
                                         ref_nones,
                                         n_haps1, n_haps2,
                                         ','.join(str(s) for s in all_haps) if len(s_star_snps) > 0 else '.'))

                                         # sum([s['arc_match'] for s in ind1_snps]),
        
        if opts.debug or local_debug: print
        if opts.debug or local_debug: print
        if opts.debug or local_debug: print
        
        pass
    return


def finish_analysis(opts):
    return
