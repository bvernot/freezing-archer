from __future__ import division
import sys, traceback, os, gzip, random
from array import array
import h5py
import argparse
from operator import itemgetter
from itertools import groupby

import time
start_time = time.time()
debug_ms = False

import locale
locale.setlocale(locale.LC_ALL, 'en_US')





def parse_arguments():

    parser = argparse.ArgumentParser(description='Manage binary representations of bed and fasta files.')
    parser.add_argument('-p', '--progress', type=int, default=50000000)
    parser.add_argument('-ref', '--reference-version', default='hg19')

    subparsers = parser.add_subparsers(dest='func')

    infile_seq_parser = argparse.ArgumentParser(add_help = False)
    infile_seq_parser.add_argument('-fa', '--fasta-files', type=argparse.FileType('r'), default=None, nargs='+', required=True)

    infile_bed_parser = argparse.ArgumentParser(add_help = False)
    infile_bed_parser.add_argument('-bed', '--bed-files', type=argparse.FileType('r'), default=None, nargs='+', required=True)
    # infile_bed_parser.add_argument('-fa', '--fasta-files', type=argparse.FileType('r'), default=None, nargs='+', required=False)
    # parser.add_argument('-vcf', '--vcf-file', required=False, type=argparse.FileType('r'), default=None)
    # parser.add_argument('-vcfz', '--gzip-vcf-file', required=False, default=None)

    bed_parser = argparse.ArgumentParser(add_help = False)
    bed_parser.add_argument('-int-chrs', '--int-chrs', default=False)
    
    compression_parser = argparse.ArgumentParser(add_help = False)
    compression_parser.add_argument('-clev', '--compression-level', type=int, default=4)
    compression_parser.add_argument('-ctype', '--compression-type', default='gzip')
    compression_parser.add_argument('-o', '--output-file', required=True, default=None)

    binary_infile_parser = argparse.ArgumentParser(add_help = False)
    binary_infile_parser.add_argument('infile')

    filetype_parser = argparse.ArgumentParser(add_help = False)
    filetype_parser.add_argument('filetype', choices=['seq', 'bed'])

    create_seq_file_parser = subparsers.add_parser('create-seq', parents = [infile_seq_parser, compression_parser])
    create_bed_file_parser = subparsers.add_parser('create-bed', parents = [infile_bed_parser, compression_parser, bed_parser])
    cat_file_parser = subparsers.add_parser('cat', parents = [binary_infile_parser])
    peek_file_parser = subparsers.add_parser('peek', parents = [binary_infile_parser])
    
    opts = parser.parse_args()
    return opts



def create_compressed_sequence_file_from_fasta(opts):

    ## create hdf5 file
    f = h5py.File(opts.output_file, 'w')
    f.attrs['version'] = myBedTools.version_str
    f.attrs['filetype'] = 'seq'
    f.attrs['ref_version'] = opts.reference_version
    
    for fa in opts.fasta_files:
    
        header = fa.readline().strip()
        if not header.startswith('>'):
            print 'INCORRECT FASTA HEADER', l
            sys.exit(-1)
            pass
        chrname = header[1:].strip().split()[0]
        print 'Extracted chromosome name "%s" from header:' % chrname
        print header
        
        ## get number of bases in chr - not a great way to do this, but.. you only have to do it once
        num_bases = 0
        for l in fa:
            if l.startswith('>'):
                print "HEADER IN MIDDLE OF FASTA FILE - myBedTools4 currently only supports single-chromosome fasta files."
                sys.exit(-1)
                pass
            num_bases += len(l)-1
            pass
        print "processing %d bases.." % num_bases
        
        ## make chr hdf5 dataset
        cdef = f.create_dataset(chrname, (num_bases,), dtype='S1', compression='gzip', compression_opts=opts.compression_level)
        
        fa.seek(0)
        x=0
        fa.readline()
        current_progress = 0
        for l in fa:
            l = l.rstrip()
            f[chrname][x:x+len(l)] = list(l)
            x += len(l)
            if opts.progress > 0 and opts.progress <= x - current_progress:
                current_progress = x
                sys.stderr.write("progress.. %s %d\n" % (header, x))
                pass
            
            pass
        
        print "finished processing fasta"
        print x, "bases out of", num_bases
        print

        pass
    
    f.close()

    return



def create_compressed_bed_file_from_bed(opts):

    ## create hdf5 file
    f = h5py.File(opts.output_file, 'w')
    f.attrs['version'] = myBedTools.version_str
    f.attrs['filetype'] = 'bed'
    f.attrs['ref_version'] = opts.reference_version

    ## make sure all chromosomes are defined..
    for chrname, num_bases in refDefinitions.chr_lens.items():
        print "Initializing chromosome %s, with %d bases" % (chrname, num_bases)
        cdef = f.create_dataset(chrname, (num_bases,), dtype='b', compression='gzip', compression_opts=opts.compression_level, fillvalue=0)
        # cdef[:] = 0
        pass

    ## find out if the chromosomes are 1 or chr1
    ## THIS DOESN'T WORK WITH STDIN
    int_chrs = False
    test_line = opts.bed_files[0].readline()
    test_chr = test_line.split()[0]
    if test_chr not in refDefinitions.chr_lens and ('chr' + test_chr) in refDefinitions.chr_lens:
        int_chrs = True
    elif test_chr not in refDefinitions.chr_lens:
        print "Can't find chromosome in hg19 reference?"
        print test_chr
        print test_line
        sys.exit(-1)
        pass
    opts.bed_files[0].seek(0)

    current_progress = 0
    for bedfile in opts.bed_files:
        for l in bedfile:
            l = l.rstrip().split()
            start = int(l[1])
            end = int(l[2])
            chr = 'chr' + l[0] if int_chrs else l[0]

            f[chr][start:end] = 1

            if opts.progress > 0 and opts.progress <= start - current_progress:
                current_progress = start
                sys.stderr.write("progress.. %s %d\n" % (chrname, start))
                pass
            pass
        pass
    
    f.close()
    
    return



def to_bed(chrom, limit=None):

    pos = 0
    count = 0

    chrname = chrom.name

    if limit != None:
        chrom = chrom[:10000]
        pass
    
    for k,g in groupby(chrom):
        l = sum(1 for x in g)
        
        if k == 0:
            pos += sum(1 for x in g)
            continue
        
        print "%s\t%d\t%d" % (chrname, pos, pos+l)
        pos += l

        count += 1
        if limit != None and count >= limit: break
        pass
    
    if limit != None and count < limit:
        # print "no elements found within 10000 bases"
        pass

    pass


def peek_compressed_file(opts):
    
    ## open hdf5 dataset
    f = h5py.File(opts.infile, 'r')
    
    if not 'version' in f.attrs:
        print "Are you sure this is the correct file?  It doesn't contain a version number.."
        sys.exit(-1)
        pass
    
    print "Version:", f.attrs['version']
    print "File contains %d chromosomes:" % len(f)

    if f.attrs['filetype'] == 'seq':
        
        for chrom in f.values():
            report_items = [chrom.name, len(chrom)] + [''.join(chrom[x:x+10]) for x in (0, len(chrom)//2, len(chrom)-100)] + [''.join(chrom[x] for x in random.sample(xrange(len(chrom)), 20))]
            print "  name:%s \t len:%d \t beginning/middle/end:%s..%s..%s \t random:%s" % tuple(report_items)
            pass
        
    elif f.attrs['filetype'] == 'bed':

        for chrom in f.values():

            report_items = [chrom.name, len(chrom)] + [''.join(str(s) for s in chrom[x:x+10]) 
                                                       for x in (0, len(chrom)//2, len(chrom)-100)] + [''.join(str(chrom[x])
                                                                                                               for x in random.sample(xrange(len(chrom)), 20))]
            print "  name:%s \t len:%d \t beginning/middle/end:%s..%s..%s \t random:%s" % tuple(report_items)

            to_bed(chrom, limit=2)

            pass
        pass


    return


class refDefinitions:

    chr_lens = {'chr1': 249250621,
                'chr2': 243199373,
                'chr3': 198022430,
                'chr4': 191154276,
                'chr5': 180915260,
                'chr6': 171115067,
                'chr7': 159138663,
                'chr8': 146364022,
                'chr9': 141213431,
                'chr10': 135534747,
                'chr11': 135006516,
                'chr12': 133851895,
                'chr13': 115169878,
                'chr14': 107349540,
                'chr15': 102531392,
                'chr16': 90354753,
                'chr17': 81195210,
                'chr18': 78077248,
                'chr19': 59128983,
                'chr20': 63025520,
                'chr21': 48129895,
                'chr22': 51304566,
                'chrX' : 155270560,
                'chrY' : 59373566,
                'chrM' : 16569}

    pass


class myBedTools:

    version_str = '0.1'

    pass



    
if __name__ == "__main__":

    opts = parse_arguments()

    if opts.func == 'create-seq':
        create_compressed_sequence_file_from_fasta(opts)
    elif opts.func == 'create-bed':
        create_compressed_bed_file_from_bed(opts)
    elif opts.func == 'cat':
        cat_compressed_file(opts)
    elif opts.func == 'peek':
        peek_compressed_file(opts)
        pass

    

    pass
