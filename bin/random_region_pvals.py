from __future__ import division
import sys
from numpy import std
import itertools
from collections import defaultdict, Counter
import arc_match_pval_tables
from time import time
from random import randint




local_debug = False

def initialize_analysis(opts):
    pass


def calc_table_match_pval(chrom, snps, opts, window_neand_pos, ind, hap1_range, hap2_range):

    # get the list of snps on each haplotype
    ind_hap1 = [s for s in snps if s['haplotypes_1'][ind] > 0 and hap1_range[0] <= s['pos'] <= hap1_range[1]]
    ind_hap2 = [s for s in snps if s['haplotypes_2'][ind] > 0 and hap2_range[0] <= s['pos'] <= hap2_range[1]]

    ret_p = [None, None]
    ret_m = [None, None]
    ret_n = [None, None]
    ret_f = [None, None]
    ret_q = [[None, None, None, None, None], [None, None, None, None, None]]

    for hapnum, ind_snps in enumerate((ind_hap1, ind_hap2)):
        
        if len(ind_snps) < 2:
            continue

        # get neand snps for this region
        neand_pos = [p for p in window_neand_pos if ind_snps[0]['pos'] <= p <= ind_snps[-1]['pos']]
        neand_sites = len(neand_pos)
        
        (test_mapped_bases_bin, test_len, \
             test_mh_sites, test_tot_sites, \
             orig_sfs, _, test_match) = arc_match_pval_tables.get_match_stats(chrom, ind_snps, neand_pos, opts, pop = 'sfs_target')

        ## HACK!  should just save sfs as an int 1-100?
        ## map test_sfs to be out of 216 (108 YRI inds in table)
        # print test_sfs, opts.num_target, int(test_sfs / opts.num_target / 2 * 216)
        test_sfs = int(orig_sfs / opts.num_target * opts.ptables_ninds_for_sfs * opts.ptables_adjust_sfs_from_target_to_ref)

        #if test_len < 10000 or test_mapped_bases_bin < 5000 | test_tot_sites < 8:
        if test_len < 10000 or test_mapped_bases_bin < 5000:
            continue
        
        eps_len = opts.len_eps
        eps_mapped_bases_bin = opts.mapped_eps
        eps_mh_sites = 1
        eps_sfs = 5
        
        factor = 1
        tot_count = 0

        if opts.table_pval_mode == 'pytables':

            current_time = time()
            match_pct = test_match / test_tot_sites
            
            query_args = {'len' : '(%f <= hap_len) & (hap_len <= %f)' % (test_len - factor * eps_len, test_len + factor * eps_len),
                          'mapped' : '(%f <= mapped_bases_bin) & (mapped_bases_bin <= %f)' % (test_mapped_bases_bin - factor * eps_mapped_bases_bin, 
                                                                                            test_mapped_bases_bin + factor * eps_mapped_bases_bin),
                          'mh' : '(%f <= mh_sites) & (mh_sites <= %f)' % (test_mh_sites - factor * eps_mh_sites, test_mh_sites + factor * eps_mh_sites),
                          'sfs' : '(%f <= sfs) & (sfs <= %f)' % (test_sfs - factor * eps_sfs, test_sfs + factor * eps_sfs)}
            
            query = ' & '.join(query_args[x] for x in opts.table_query)

            if test_match > 0:
                a = [ x[:] for x in opts.match_pval_table.root.table.where(query) ]
            else:
                ## no need to actually query for similar haps, if we know that the pval is going to be 1
                a = []
                pass

            # len > 10000 & mapped_bases_bin > 5000 & tot_sites >= 8
            
            # names = ['count', 'mapped_bases_bin', 'len', 'mh_sites', 'tot_sites', 'sfs', 'std_dev', 'match'])
            # 1 10000 11111 20 28 141 62 14.0
            # 1 10000 11111 20 32 106 82 10.0
            # 6 10000 11111 21 36 170 51 15.0
            
            tot_count = sum(x[0] for x in a)
            
            match_pct_by_rows = [x[0] if x[7]/x[4] >= match_pct else 0 for x in a]
            pval = (sum(match_pct_by_rows) + 1) / (tot_count + 1)
            print 'PYTABLES', match_pct, pval, len(a), tot_count, time()-current_time, query

        else:
            
            query_args = {'len' : '((@test_len-@factor*@eps_len) <= len <= (@test_len+@factor*@eps_len))',
                          'mapped' : '((@test_mapped_bases_bin-@factor*@eps_mapped_bases_bin) <= mapped_bases_bin ' + \
                              '<= (@test_mapped_bases_bin+@factor*@eps_mapped_bases_bin))',
                          'mh' : '((@test_mh_sites-@factor*@eps_mh_sites) <= mh_sites <= (@test_mh_sites+@factor*@eps_mh_sites))',
                          'sfs' : '((@test_sfs-@factor*@eps_sfs) <= sfs <= (@test_sfs+@factor*@eps_sfs))'}
            
            query = ' & '.join(query_args[x] for x in opts.table_query)
            
            # print 'QEURY', query, test_len, test_mapped_bases_bin, test_mh_sites, test_sfs
            
            tot_count = 0

            # while tot_count < 1000:
                
            current_time = time()
            t = opts.match_pval_table[0].query(query, engine='numexpr')
            fail_count = sum(t.query('match / tot_sites >= @test_match / @test_tot_sites').eval('count'))
            tot_count  = sum(t['count'])

            # print 'QUERY', tot_count, factor, fail_count
            # print test_len-factor*eps_len, test_len+factor*eps_len
            # print test_mapped_bases_bin-factor*eps_mapped_bases_bin, test_mapped_bases_bin+factor*eps_mapped_bases_bin
            # print test_mh_sites-factor*eps_mh_sites, test_mh_sites+factor*eps_mh_sites
            # print test_sfs-factor*eps_sfs, test_sfs+factor*eps_sfs

            # factor *= 1.1

            pval = (fail_count+1)/(tot_count+1)
            match_pct = test_match / test_tot_sites
            
            print 'PANDAS', match_pct, pval, tot_count, time()-current_time, query

            pass

        ret_p[hapnum] = pval
        ret_m[hapnum] = match_pct
        ret_n[hapnum] = tot_count
        ret_f[hapnum] = factor / 1.1
        ret_q[hapnum] = (test_len, test_mapped_bases_bin, test_mh_sites, test_sfs, orig_sfs)
        pass
        
    return (ret_p, ret_m, ret_n, ret_f, ret_q)


def run_window_analysis(chrom, winstart, winend, snps, opts):

    if opts.debug: print 'starting window', chrom, winstart, winend, 'NUM SNPS:', len(snps)
    if len(snps) <= 2: return
    
    window_neand_pos = opts.archaic_vcf.get_derived_sites(chrom, winstart, winend)


    ## JUST LOOP OVER TARGET INDS
    for ind in opts.target_indices:

        # only consider 1/1000 possible pvals?
        if randint(0,99) != 99: continue

        ind1_snps = [s for s in snps if s['genotypes'][ind] > 0 and s['target'] and not s['reference']]
        ind1_snps_all = [s for s in snps if s['genotypes'][ind] > 0]

        # used as a count of the number of snvs in a region, for matching to the null model
        ind1_or_ref_snps = [s for s in snps if s['genotypes'][ind] > 0 or s['reference']]

        ind1_pos = [s['pos'] for s in ind1_snps]

        # ind1_hap1 and ind1_hap2 are just masks for which snps are on which haplotypes - they're the same length as ind1_snps
        ind1_hap1 = [s['haplotypes_1'][ind] > 0 for s in ind1_snps]
        ind1_hap2 = [s['haplotypes_2'][ind] > 0 for s in ind1_snps]


        ## get S* haplotype range
        hap1_pos = [s['pos'] for s in ind1_snps if s['haplotypes_1'][ind] > 0]
        #hap1_pos = [ind1_snps[i]['pos'] for i in s_star_snps if ind1_hap1[i]]
        hap1_range = (min(hap1_pos), max(hap1_pos)) if len(hap1_pos) > 1 else (0,0)
        hap2_pos = [s['pos'] for s in ind1_snps if s['haplotypes_2'][ind] > 0]
        #hap2_pos = [ind1_snps[i]['pos'] for i in s_star_snps if ind1_hap2[i]]
        hap2_range = (min(hap2_pos), max(hap2_pos)) if len(hap2_pos) > 1 else (0,0)

        

        if opts.match_pval_table != None:
            (match_table_pvals, match_table_pcts, \
                 match_table_Ns, match_table_factors, match_table_query) = calc_table_match_pval(chrom, snps, opts, window_neand_pos, ind, \
                                                                                                     hap1_range, hap2_range)
        else:
            match_table_pvals = (None, None)
            match_table_pcts = (None, None)
            match_table_Ns = (None, None)
            match_table_factors = (None, None)
            match_table_query = [[None, None, None, None, None], [None, None, None, None, None]]
            pass


        ms_fields = []
        ms_fields_labels = []
        if opts.vcf_is_ms_file:
            
            # count sites that are intr / S* / both

            # intr and S*: sites that are a) in the sstar hap range, b) introgressed on that hap, and c) on that hap at all (this may be redundant!)
            ms_intr_ss_sites_hap1 = len([s for s in ind1_snps_all if (hap1_range[0] <= s['pos'] <= hap1_range[1]) and s['haplotypes_1_intr'][ind] and s['haplotypes_1'][ind] > 0])
            ms_intr_ss_sites_hap2 = len([s for s in ind1_snps_all if (hap2_range[0] <= s['pos'] <= hap2_range[1]) and s['haplotypes_2_intr'][ind] and s['haplotypes_2'][ind] > 0])

            # S*: sites that are a) in the sstar hap range, and b) on that hap at all
            ms_ss_sites_hap1      = len([s for s in ind1_snps_all if (hap1_range[0] <= s['pos'] <= hap1_range[1]) and s['haplotypes_1'][ind] > 0])
            ms_ss_sites_hap2      = len([s for s in ind1_snps_all if (hap2_range[0] <= s['pos'] <= hap2_range[1]) and s['haplotypes_2'][ind] > 0])

            # intr: sites that are a) introgressed on that hap, and c) on that hap at all (this may be redundant!)
            ms_intr_sites_hap1      = len([s for s in ind1_snps_all if s['haplotypes_1_intr'][ind] and s['haplotypes_1'][ind] > 0])
            ms_intr_sites_hap2      = len([s for s in ind1_snps_all if s['haplotypes_2_intr'][ind] and s['haplotypes_2'][ind] > 0])

            # N: sites that are a) on that hap at all
            ms_N_sites_hap1      = len([s for s in ind1_snps_all if s['haplotypes_1'][ind] > 0])
            ms_N_sites_hap2      = len([s for s in ind1_snps_all if s['haplotypes_2'][ind] > 0])

            # print "ind1_snps_all hap1", [s['pos'] for s in ind1_snps_all if s['haplotypes_1'][ind] > 0], 'intr:', [s['haplotypes_1_intr'][ind] for s in ind1_snps_all if s['haplotypes_1'][ind] > 0]
            # print "ind1_snps_all hap2", [s['pos'] for s in ind1_snps_all if s['haplotypes_2'][ind] > 0], 'intr:', [s['haplotypes_2_intr'][ind] for s in ind1_snps_all if s['haplotypes_2'][ind] > 0]
            
            ms_fields += [ms_intr_sites_hap1,
                          ms_intr_sites_hap2,
                          ms_ss_sites_hap1,
                          ms_ss_sites_hap2,
                          ms_intr_ss_sites_hap1,
                          ms_intr_ss_sites_hap2,
                          ms_N_sites_hap1,
                          ms_N_sites_hap2]

            ms_fields_labels += ['ms_intr_sites_hap1',
                                 'ms_intr_sites_hap2',
                                 'ms_ss_sites_hap1',
                                 'ms_ss_sites_hap2',
                                 'ms_intr_ss_sites_hap1',
                                 'ms_intr_ss_sites_hap2',
                                 'ms_N_sites_hap1',
                                 'ms_N_sites_hap2']# ,
                                 # 'ms_intr_bases_hap1',
                                 # 'ms_intr_bases_hap2',
                                 # 'ms_ss_bases_hap1',
                                 # 'ms_ss_bases_hap2',
                                 # 'ms_intr_ss_bases_hap1',
                                 # 'ms_intr_ss_bases_hap2']
            
            pass

        if opts.first_line:
            print '\t'.join(['chrom', 'winstart', 'winend', 'n_snps', 'n_ind_snps', 'n_region_ind_snps',
                             'ind_id',
                             'pop',
                             'hap_1_window_pval_table', 'hap_2_window_pval_table',
                             'hap_1_window_match_pct_table', 'hap_2_window_match_pct_table',
                             'hap_1_window_match_N_table', 'hap_2_window_match_N_table',
                             'hap_1_window_match_f_table', 'hap_2_window_match_f_table',
                             'hap_1_window_match_len_table', 'hap_2_window_match_len_table',
                             'hap_1_window_match_mapped_table', 'hap_2_window_match_mapped_table',
                             'hap_1_window_match_mh_table', 'hap_2_window_match_mh_table',
                             'hap_1_window_match_sfs_table', 'hap_2_window_match_sfs_table',
                             'hap_1_window_match_orig_sfs_table', 'hap_2_window_match_orig_sfs_table',
                             'hap_1_start', 'hap_1_end', 
                             'hap_2_start', 'hap_2_end'] + ms_fields_labels +
                            opts.tag_ids)

            opts.first_line = False
            pass

        print '\t'.join(str(s) for s in [chrom, winstart, winend, len(snps), len(ind1_snps), len(ind1_or_ref_snps),
                                         opts.get_id_from_sample_index(ind),
                                         opts.get_pop_from_sample_index(ind),
                                         match_table_pvals[0], match_table_pvals[1],
                                         match_table_pcts[0], match_table_pcts[1],
                                         match_table_Ns[0], match_table_Ns[1],
                                         match_table_factors[0], match_table_factors[1],
                                         match_table_query[0][0], match_table_query[1][0],
                                         match_table_query[0][1], match_table_query[1][1],
                                         match_table_query[0][2], match_table_query[1][2],
                                         match_table_query[0][3], match_table_query[1][3],
                                         match_table_query[0][4], match_table_query[1][4],
                                         hap1_range[0], hap1_range[1],
                                         hap2_range[0], hap2_range[1]] +
                        ms_fields +
                        opts.tags)

                                         # sum([s['arc_match'] for s in ind1_snps]),
        
        if opts.debug or local_debug: print
        if opts.debug or local_debug: print
        if opts.debug or local_debug: print
        
        pass
    return


def finish_analysis(opts):
    return
