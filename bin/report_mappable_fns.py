from __future__ import division
import sys
from numpy import std
import itertools
from collections import defaultdict, Counter


def initialize_analysis(opts):
    pass


def run_window_analysis(chrom, winstart, winend, snps, opts):
    # print 'starting window', chrom, winstart, winend, 'NUM SNPS:', len(snps)

    full_chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + chrom
    my_mapped_bases = opts.regions.amount_in_region(full_chrom, winstart, winend)

    print '\t'.join(str(s) for s in [chrom, winstart, winend, my_mapped_bases, len(snps)])
        
    return


def finish_analysis(opts):
    return
