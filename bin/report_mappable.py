from __future__ import division
import sys
from numpy import std
import itertools
from collections import defaultdict, Counter


def initialize_analysis(opts):
    pass


def run_window_analysis(chrom, winstart, winend, snps, opts):
    # print 'starting window', chrom, winstart, winend, 'NUM SNPS:', len(snps)

    for s in snps:
        #print s

        if sum([s['genotypes'][i] for i in opts.target_indices]) != s['sfs_target']:
            print "error in sfs_target?"
            print s
            sys.exit(-1)
            pass

        if sum([s['genotypes'][i] for i in opts.reference_indices]) != s['sfs_reference']:
            print "error in sfs_reference?"
            print s
            print sum([s['genotypes'][i] for i in opts.reference_indices]), '!=', s['sfs_reference']
            print 'sum(', [s['genotypes'][i] for i in opts.reference_indices], ')', '!=', s['sfs_reference']
            print opts.reference_indices
            print opts.num_target, opts.num_reference
            sys.exit(-1)
            pass

        print '\t'.join(str(x) for x in [s['chrom'], s['pos'], \
                                             s['sfs_target'], s['sfs_reference']] + \
                            ['%d|%d' % (s['haplotypes_1'][i], s['haplotypes_2'][i]) for i in opts.target_indices + opts.reference_indices] + \
                            [opts.archaic_vcf.get_derived_count(s['chrom'], s['pos']), s['arc_match']])
        pass

    print "archaic derived sites"
    for p, d in opts.archaic_vcf.get_derived_sites_with_der_count(chrom, winstart, winend):
        print '\t'.join(str(x) for x in [chrom, p, d])
        pass
        
    return


def finish_analysis(opts):
    return
