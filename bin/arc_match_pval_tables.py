from __future__ import division
from numpy import std
from collections import defaultdict, Counter
import sys

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


def get_match_stats(chrom, ind_snps, neand_pos, opts, pop = 'sfs_reference'):
    
    # get the total number of sites on this haplotype
    ind_pos = set([s['pos'] for s in ind_snps])
    mh_sites = len(ind_pos)

    # get the number of mapped bases in this region
    mapped_bases = opts.regions.amount_in_region(chrom, min(ind_pos), max(ind_pos)) if opts.regions != None else max(ind_pos) - min(ind_pos)
    mapped_bases_bin = mapped_bases // 1000 * 1000
    
    # get the total number of sites we are "worrying" about - i.e., the number of places where N and MH could diverge
    ind_pos2 = ind_pos.union(set(neand_pos))
    tot_sites = len(ind_pos2)

    # count the number of matches to neanderthal
    match_sites = sum([s['arc_match'] for s in ind_snps])
    match_dist  = sum([opts.archaic_vcf.get_derived_count(s['chrom'], s['pos'])/2 for s in ind_snps if s['arc_match']])
    
    # get the average freq
    sfs = int(round(sum([s[pop] for s in ind_snps]) / mh_sites))
    std_dev = int(round(std([s[pop] for s in ind_snps])))
    
    return (mapped_bases_bin, max(ind_pos) - min(ind_pos), \
                mh_sites, tot_sites, \
                sfs, std_dev, match_dist)

def run_window_analysis(chrom, winstart, winend, snps, opts):

    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom

    winends = list(reversed(range(winstart+10000, winend, 1000)))

    #print winends

    for winend in winends:

        ## first remove any snps past the new winend
        snps = [s for s in snps if s['pos'] < winend and s['reference']]
        #print len(snps)
        if len(snps) == 0: continue

        snp_pos = set([s['pos'] for s in snps])

        if max(snp_pos) < winend-1000:
            #print 'SKIPPING WINEND, LAST SNP DOES NOT END IN THIS WINDOW', winend, winend-1000, max(snp_pos), snp_pos
            continue
        #print 'KEEPING  WINEND, LAST SNP DOES     END IN THIS WINDOW', winend, winend-1000, max(snp_pos), snp_pos
        

        # get neand snps for this region
        neand_pos = opts.archaic_vcf.get_derived_sites(chrom, winstart, winend)
        neand_sites = len(neand_pos)

        # JUST LOOKING AT REFERENCE INDIVIDUALS.  I.E., YRI FROM S* CALCULATION
        for ind in opts.reference_indices:

            # get the list of genotypes for each individual (ignore target and reference?)  ## and s['target'] and not s['reference']
            # ind_snps = [s for s in snps if s['genotypes'][ind] > 0]
            ind_hap1 = [s for s in snps if s['haplotypes_1'][ind] > 0]
            ind_hap2 = [s for s in snps if s['haplotypes_2'][ind] > 0]

            for hapnum, ind_snps in enumerate((ind_hap1, ind_hap2)):

                if len(ind_snps) < 1: continue

                # make sure our last snp is in the window we're considering
                if ind_snps[-1]['pos'] < winend-1000:
                    # print 'SKIPPING HAP, DOES NOT END IN THIS WINDOW', winend, winend-1000, max(ind_pos), ind_pos
                    continue
                # print 'KEEPING  HAP, DOES     END IN THIS WINDOW', winend, winend-1000, max(ind_pos), ind_pos
                
                (mapped_bases_bin, hap_length, \
                     mh_sites, tot_sites, \
                     sfs, std_dev, match_dist) = get_match_stats(full_chrom, ind_snps, neand_pos, opts)

                print 'STATS', chrom, winstart, winend, \
                    opts.get_id_from_sample_index(ind), opts.get_pop_from_sample_index(ind), hapnum, \
                    mapped_bases_bin, hap_length, \
                    mh_sites, neand_sites, tot_sites, sfs, std_dev, match_dist

                opts.arc_match_table[mapped_bases_bin][tot_sites][sfs][std_dev].update([match_dist])

                if opts.debug or local_debug: print
                pass

            pass
        pass
    return


def finish_analysis(opts):
    # # print opts.arc_match_table
    # import cPickle

    #print opts.arc_match_table
    print_lim = 10
    print_lim = sys.maxint

    for i in opts.arc_match_table.keys()[:print_lim]:
        #print i
        for j in opts.arc_match_table[i].keys()[:print_lim]:
            #print i,j
            for k in opts.arc_match_table[i][j].keys()[:print_lim]:
                #print i,j,k
                for l in opts.arc_match_table[i][j][k].keys()[:print_lim]:
                    for m in opts.arc_match_table[i][j][k][l].keys()[:print_lim]:
                        print 'TBL', i,j,k,l,m,opts.arc_match_table[i][j][k][l][m]
                        pass
                    pass
                pass
            pass
        pass

    # cPickle.dump(opts.arc_match_table, opts.pickle_file_for_match_pvals, -1)
    # opts.pickle_file_for_match_pvals.close()
    return
