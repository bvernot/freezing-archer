from __future__ import division
from numpy import std
from collections import defaultdict, Counter

local_debug = False

def mydd():
    return defaultdict( Counter )
def mydd2():
    return defaultdict( mydd )
def mydd3():
    return defaultdict( mydd2 )

def initialize_analysis(opts):
    
    print "initializing tables..."
    #opts.arc_match_table = defaultdict(lambda : defaultdict(lambda : defaultdict( Counter )))
    opts.arc_match_table = defaultdict(mydd3)

    # print opts.arc_match_table[0][0][0]
                
    print "... FINISHED initializing tables"
        
    pass


def run_window_analysis(chrom, winstart, winend, snps, opts):

    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    mapped_bases = opts.regions.amount_in_region(full_chrom, winstart, winend)
    mapped_bases_bin = mapped_bases // 1000 * 1000

    # get neand snps for this region
    neand_pos = [p for p in range(winstart, winend) \
                     if opts.neand_vcf.has_var_at_site(chrom, p) and \
                     not opts.neand_vcf.has_ref(chrom, p) and \
                     (opts.regions == None or opts.regions.in_region_one_based(full_chrom, p))]
    # print 'neand positions', chrom, winstart, winend, \
    #    mapped_bases, len(snps), len(opts.neand_vcf.vcf[chrom]), len(neand_pos), neand_pos, [p for p in range(winstart, winend) if opts.neand_vcf.has_var_at_site(chrom, p)]

    for ind in range(len(opts.updated_ind_ids)):

        # get the list of genotypes for each individual (ignore target and reference?)  ## and s['target'] and not s['reference']
        ind_snps = [s for s in snps if s['genotypes'][ind] > 0]
        if len(ind_snps) == 0: continue

        # get the total number of sites we are "worrying" about - i.e., the number of places where N and MH could diverge
        ind_pos = set([s['pos'] for s in ind_snps])
        ind_pos2 = ind_pos.union(set(neand_pos))
        n_sites = len(ind_pos2)
        #print len(ind_pos), len(neand_pos), len(ind_pos2)
        #print neand_pos
        #print ind_pos
        #print ind_pos2
        #print

        # count the number of matches to neanderthal
        n_match = sum([s['arc_match'] for s in ind_snps])
        # for s in ind_snps:
        #     print 'nmatchdebug', chrom, winstart, winend, opts.updated_ind_ids[ind], n_sites, len(ind_snps), n_match, s['pos'], s['arc_match']
        #     pass
        
        # get the average freq
        sfs = int(round(sum([s['sfs_target'] for s in ind_snps]) / n_sites))

        # get the std freq
        std_dev = int(round(std([s['sfs_target'] for s in ind_snps])))

        if opts.debug or local_debug: print chrom, winstart, winend, n_match, n_sites, sfs, std_dev
        if opts.debug or local_debug: print [s['sfs_target'] for s in ind_snps]
        if opts.debug or local_debug: print [s['genotypes'][ind] for s in ind_snps]

        # store in the table
        # print chrom, winstart, winend, ind, opts.arc_match_table[mapped_bases_bin][n_sites][sfs][std_dev], n_match
        opts.arc_match_table[mapped_bases_bin][n_sites][sfs][std_dev].update([n_match])

        if opts.debug or local_debug: print
        
        pass
    return


def finish_analysis(opts):
    # print opts.arc_match_table
    import cPickle

    print_lim = 2

    for i in opts.arc_match_table.keys()[:print_lim]:
        print i
        for j in opts.arc_match_table[i].keys()[:print_lim]:
            print i,j
            for k in opts.arc_match_table[i][j].keys()[:print_lim]:
                print i,j,k
                for l in opts.arc_match_table[i][j][k].keys()[:print_lim]:
                    print i,j,k,l,opts.arc_match_table[i][j][k][l]
                    pass
                pass
            pass
        pass
    cPickle.dump(opts.arc_match_table, opts.pickle_file_for_match_pvals, -1)
    opts.pickle_file_for_match_pvals.close()
    return
