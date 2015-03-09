from __future__ import division
import sys, itertools
import argparse
from read_ms import process_ms_block_to_genotypes


import time
start_time = time.time()

import locale
locale.setlocale(locale.LC_ALL, 'en_US')



parser = argparse.ArgumentParser(description='Convert an ms file to a VCF file - in a very basic way.')
parser.add_argument('msfile', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--output-file', type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', '--debug', action='store_true')
parser.add_argument('-winlen', '--window-length', type=int, default=50000)
parser.add_argument('-winstep', '--window-step', type=int, default=20000)

opts = parser.parse_args()

# ms 42 1 -s 10 
# 44126 40565 42561

# //
# segsites: 10
# positions: 0.1717 0.2230 0.2277 0.4523 0.4598 0.5201 0.7094 0.8533 0.9100 0.9894 
# 0010000001
# 0110000001
# 1000000000
# 1000000000
# 0010000001

# header
header = opts.msfile.readline().strip().split()

# parse params from header
(_, nsam, howmany) = header[:3]
nsam = int(nsam)
num_inds = nsam // 2
#setattr(opts, 'num_inds', int(nsam) // 2)

if nsam != 2 * num_inds:
    sys.stderr.write('Not an even number of chromosomes - dropping last haplotype.\n')
    pass

## print header
print '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + ['i%d' % d for d in xrange(1,num_inds+1)])

# seed params
seed_params = opts.msfile.readline()

line = opts.msfile.readline()
while not line.startswith('//'):
    if opts.debug: print line
    line = opts.msfile.readline()
    pass

snps = process_ms_block_to_genotypes(opts.msfile, num_inds, opts.window_length)

if opts.debug: 
    snps = list( snps )
    print snps
    pass

blocknum = 1
for pos, gts in snps:
    print '\t'.join(str(s) for s in ['ms%d' % blocknum, pos, '.', 'A', 'G', '.', '.', '.', 'GT'] + ['%d|%d' % (gt[0], gt[1]) for gt in gts])
    pass


