from __future__ import division
import sys, traceback, os, gzip
from array import array
import h5py
import argparse
from operator import itemgetter

import time
start_time = time.time()
debug_ms = False

import locale
locale.setlocale(locale.LC_ALL, 'en_US')



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Manage binary representations of bed and fasta files.')

    # parser.add_subparsers
    
    # parser.add_argument('-vcf', '--vcf-file', required=False, type=argparse.FileType('r'), default=None)
    parser.add_argument('-fa', '--fasta-file', required=False, type=argparse.FileType('r'), default=None)
    # parser.add_argument('-vcfz', '--gzip-vcf-file', required=False, default=None)
    # parser.add_argument('-indf', '--ind-pop-file', required=True, type=argparse.FileType('r'))
    
    parser.add_argument('-o', '--output-file', required=False, default=None)
    parser.add_argument('-p', '--progress', type=int, default=5000000)
    parser.add_argument('-clev', '--compression-level', type=int, default=4)

    opts = parser.parse_args()

    header = opts.fasta_file.readline().strip()
    if not header.startswith('>'):
        print 'INCORRECT FASTA HEADER', l
        sys.exit(-1)
        pass

    ## get number of bases in chr - not a great way to do this, but.. you only have to do it once
    num_bases = 0
    for l in opts.fasta_file:
        num_bases += len(l)-1
        pass
    print "processing %d bases.." % num_bases

    ## make chr hdf5 dataset
    f = h5py.File(opts.output_file, 'w')
    cdef = f.create_dataset('c', (num_bases,), dtype='S1', compression='gzip', compression_opts=opts.compression_level)

    opts.fasta_file.seek(0);
    x=0
    opts.fasta_file.readline()
    current_progress = 0
    for l in opts.fasta_file:
        l = l.rstrip()
        f['c'][x:x+len(l)] = list(l)
        x += len(l)
        if opts.progress > 0 and opts.progress <= x - current_progress:
            current_progress = x
            sys.stdout.write("progress.. %s %d\n" % (header, x))
            pass

        pass
    
    print "finished processing fasta"
    print x, "bases out of", num_bases
    

    pass



class myBedTools:
    
    pass



    
