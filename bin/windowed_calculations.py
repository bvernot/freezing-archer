from __future__ import division
import sys, os, random, itertools, shutil, cStringIO
import argparse
from read_vcf import vcf_to_genotypes_windowed
from read_ms import ms_to_genotypes_windowed
import gzip
from custom_argparse import *
import cPickle

#from numpy import array, arange
# sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/mapped_snps/')
# sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/metrics/')
# sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/filter_files/')
# sys.path.append('/net/akey/vol1/home/bvernot/archaic_exome/experiments/fdr_simulated_basic/latest/')
# #from region_stats import region_type_stats
# from myBedTools3 import myBedTools
# from argparse_bedtools import *
# from BaseLookup import BaseLookup
# import fileinput
# from operator import itemgetter
# import argparse
# #from numpy.random import binomial
# #import tables
# #import re
# #import sqlite3
# from bitarray import bitarray
# from collections import Counter, defaultdict
# import get_pct_arc_per_ind_from_ms_file_new as pct_arc
# from mydefaultdict import mydefaultdict
# #from test_parse_tree import read_tree4, find_node, tree_to_str, get_terminals, get_path_to_root, get_dist_btwn_nodes

import time
start_time = time.time()
debug_ms = False

import locale
locale.setlocale(locale.LC_ALL, 'en_US')



parser = argparse.ArgumentParser(description='Calculate S*')

parser.add_argument('-vcf', '--vcf-file', required=False, type=argparse.FileType('r'), default=None)
parser.add_argument('-vcfz', '--gzip-vcf-file', required=False, default=None)
parser.add_argument('-indf', '--ind-pop-file', required=False, type=argparse.FileType('r'), default=None)
parser.add_argument('-ptable', '--match-pval-table', type=argparse.FileType('rb'), default=None, nargs='+')
parser.add_argument('-pfile', '--pickle-file-for-match-pvals', type=argparse.FileType('wb'), default=None)
parser.add_argument('-winlen', '--window-length', type=int, default=50000)
parser.add_argument('-winstep', '--window-step', type=int, default=20000)
parser.add_argument('-l', '--limit-wins', type=int, default=0)
parser.add_argument('-p', '--progress', type=int, default=0)
parser.add_argument('-ref-pops', '--reference-populations', default=[], nargs='+')
parser.add_argument('-ref-inds', '--reference-individuals', default=[], nargs='+')
parser.add_argument('-target-pops', '--target-populations', default=[], nargs='+')
parser.add_argument('-target-inds', '--target-individuals', default=[], nargs='+')
parser.add_argument('-exclude-pops', '--exclude-populations', default=[], nargs='+', 
                    help='Exclude snps where these individuals have a nonref allele (or derived, if an ancestral genome is given)')
parser.add_argument('-exclude-inds', '--exclude-individuals', default=[], nargs='+',
                    help='Exclude snps where these individuals have a nonref allele (or derived, if an ancestral genome is given)')
parser.add_argument('-tag-ids', '--tag-ids', default=[], nargs='+')
parser.add_argument('-tags', '--tags', default=[], nargs='+')
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-ms', '--vcf-is-ms-file', action='store_true')
parser.add_argument('-mspops', '--ms-pop-sizes', default=None, nargs='+', type=int, help='This is identical to the -I argument for ms. WRT target and reference populations, numbering starts from 0.')
parser.add_argument('-msinds', '--ms-num-diploid-inds', default=None, type=int, help='The number of diploid individuals considered. This is important because we sometimes simulate a single archaic chromosome.')
# parser.add_argument('-msarc', '--ms-archaic-chromosomes', default=None, nargs='+', type=int, help='The archaic chromosomes, if simulated.')
parser.add_argument('-msarc', '--ms-archaic-populations', default=[], nargs='+', type=int, help='The archaic populations, if simulated.')
parser.add_argument('-illumina-chrom', '--vcf-has-illumina-chrnums', action='store_true')
parser.add_argument('-archaic-vcf', '--archaic-vcf', action = VCFFileAction, required = False, nargs='+', help = 'VCF file listing archaic sites', default=None)
parser.add_argument('-ancbsg', '--ancestral-bsg', action = BinarySeqFileAction, required = False, help = 'BSG file listing ancestral sites (CAnc or just chimp)')

parser.add_argument('-r', '--regions', action = BinaryBedFileAction, required=False, default=None, 
                    help = 'A bbg file that specifies which regions to consider.  Only snps in this region are loaded.')
parser.add_argument('-ir', '--intersect-region', action = IntersectBinaryBedFilesAction, nargs='+', required=False, default=None, 
                    help = 'A bbg file that is intersected with the --regions file to produce a new set of regions to consider.')
parser.add_argument('-x', '--exclude-region', nargs='+', action = MergeBinaryBedFilesAction, required=False, default=None, 
                    help = 'bbg file(s) that specify which regions should be excluded from the analysis.  If more than one file is given, the files are merged.')
parser.add_argument('-o', '--output-file', type=argparse.FileType('w'), required = False, default = sys.stdout, help = 'output file')

analysis_group = parser.add_mutually_exclusive_group(required=True)
analysis_group.add_argument('-s-star', '--s-star', action='store_true')
analysis_group.add_argument('-test-fns', '--test-fns', action='store_true')
analysis_group.add_argument('-match-table', '--make-arc-match-pval-tables', action='store_true')
analysis_group.add_argument('-d-stats', '--d-statistics', action='store_true')



opts = parser.parse_args()
setattr(opts, 'first_line', True)

# set stdout to the -o file, if applicable
sys.stdout = opts.output_file

## select the appropriate set of functions
## this is.. probably poor form (runtime import selection), but it is easy!

if opts.s_star:
    from s_star_fns \
        import initialize_analysis, run_window_analysis, finish_analysis
elif opts.make_arc_match_pval_tables:
    from arc_match_pval_tables \
        import initialize_analysis, run_window_analysis, finish_analysis
elif opts.d_statistics:
    from d_stat_fns \
        import initialize_analysis, run_window_analysis, finish_analysis
elif opts.test_fns:
    from test_fns \
        import initialize_analysis, run_window_analysis, finish_analysis
    pass

## THIS IS SUCH A HACK - ONLY USING ONE ARCHAIC VCF AT THIS POINT, SO IF THERE'S MORE THAN ONE...... JUST REMOVE THEM
if opts.archaic_vcf != None:
    opts.archaic_vcf = opts.archaic_vcf[0]
    pass

munge_regions(opts)


if opts.vcf_file == None and opts.gzip_vcf_file == None:
    print "Require at least one vcf file option."
    sys.exit(-1)
    pass
elif opts.vcf_file != None and opts.gzip_vcf_file != None:
    print "Require exactly one vcf file option."
    sys.exit(-1)
    pass

if not opts.vcf_is_ms_file and opts.ind_pop_file == None:
    print "VCF file requires --ind-pop-file."
    sys.exit(-1)
    
elif opts.vcf_is_ms_file and opts.ms_pop_sizes == None:
    print "ms file requires --ms-pop-sizes."
    sys.exit(-1)

elif opts.vcf_is_ms_file and opts.ms_num_diploid_inds == None:
    print "ms file requires --ms-num-diploid-inds."
    sys.exit(-1)
    
# elif opts.vcf_is_ms_file and opts.ms_archaic_populations == None:
#     print "ms file requires --ms-archaic-populations."
#     sys.exit(-1)
#     pass
    

if opts.gzip_vcf_file != None:
    opts.vcf_file = gzip.open(opts.gzip_vcf_file)
    pass
            
if opts.vcf_is_ms_file:
    read_genotype_fn = ms_to_genotypes_windowed
else:
    read_genotype_fn = vcf_to_genotypes_windowed
    pass

if opts.match_pval_table != None:
    t = []
    for f in opts.match_pval_table:
        sys.stderr.write('loading pval table: %s\n' % f)
        t.append(cPickle.load(f))
        f.close()
        pass
    opts.match_pval_table = t
    sys.stderr.write('...finished loading pval tables\n')
    pass



initialize_analysis(opts)

nwins = 0
for chrom, winstart, winend, snps in read_genotype_fn(opts.vcf_file, opts.window_length, opts.window_step, opts.ind_pop_file, opts):

    nwins += 1
    local_debug = False

    if opts.debug or local_debug: print
    if opts.debug or local_debug: print winstart, winend, len(snps)
    if opts.debug or local_debug: print opts.target_indices

    if opts.debug or local_debug:
        for snp in snps:
            print 'snp', snp['pos'], snp
            pass
        pass

    run_window_analysis(chrom, winstart, winend, snps, opts)

    if opts.limit_wins > 0 and nwins >= opts.limit_wins: break
    if opts.progress > 0 and nwins % opts.progress == 0: sys.stderr.write("progress.. %s %d %d\n" % (chrom, winstart, winend))


    pass



finish_analysis(opts)
