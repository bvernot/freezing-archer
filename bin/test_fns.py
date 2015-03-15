from __future__ import division
import sys
from numpy import std
import itertools
from collections import defaultdict, Counter


def initialize_analysis(opts):
    pass


def run_window_analysis(chrom, winstart, winend, snps, opts):
    print 'starting window', chrom, winstart, winend, 'NUM SNPS:', len(snps)

    for s in snps:
        print s

    return


def finish_analysis(opts):
    return
