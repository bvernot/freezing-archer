from __future__ import division
import sys, math
sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/mapped_snps/')
from BaseLookup import BaseLookup



if len(sys.argv) < 3:
    print 'useage:'
    print '%s window_size window_step [chr start end]' % sys.argv[0]
    sys.exit(-1)
    pass

window_size = int(sys.argv[1])
window_step = int(sys.argv[2])

if len(sys.argv) > 3:
    region = True
    chr_list = [sys.argv[3]]
    start_reg = int(sys.argv[4])
    end_reg_fn = lambda c : int(sys.argv[5])
else:
    region = False
    chr_list = BaseLookup.chrs
    start_reg = 0
    end_reg_fn = lambda c : BaseLookup.chr_lens[c]
    pass

for chr in chr_list:
    start = start_reg

    if chr == 'chrX':

        while start < BaseLookup.chrX_par1 and start < end_reg_fn(chr):
            print "%s\t%d\t%d" % (chr, start, min(start + window_size, BaseLookup.chrX_par1, end_reg_fn(chr)))
            start += window_step
            pass

        start = BaseLookup.chrX_par1
        while start < BaseLookup.chrX_par2 and start < end_reg_fn(chr):
            print "%s\t%d\t%d" % (chr, start, min(start + window_size, BaseLookup.chrX_par2, end_reg_fn(chr)))
            start += window_step
            pass

        start = BaseLookup.chrX_par2
        while start < end_reg_fn(chr):
            print "%s\t%d\t%d" % (chr, start, min(start + window_size, end_reg_fn(chr)))
            start += window_step
            pass
        
        pass
    else:
        while start < end_reg_fn(chr):
            print "%s\t%d\t%d" % (chr, start, min(start + window_size, end_reg_fn(chr)))
            start += window_step
            pass
        pass
    pass

