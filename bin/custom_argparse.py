import argparse, sys
from collections import defaultdict
import locale
import gzip
import gc

from myBedTools3 import myBedTools

## in original mybedtools, I allowed the user to specify the ref version with namespace.ref_version
## here I've replaced all of those with:
local_ref_version = 'b37'

class BinarySeqFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.binaryseqfilegenome, initialize = False, ref_version = local_ref_version)
        bt.read_as_regions(filename, file_type=myBedTools.binaryseqfilegenome)
        setattr(namespace, self.dest, bt)
        pass
    pass


class missingdict(defaultdict):
    # missingdict allows us to not add something to the defaultdict just because we query it
    def __missing__(self, key):
        return self.default_factory()
    pass



class BinaryBedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = local_ref_version)
        bt.read_as_regions(filename)
        setattr(namespace, self.dest, bt)
        pass
    pass

class MultBinaryBedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        def fn_outer(i):
            bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = local_ref_version)
            bt.read_as_regions(i)
            return bt
        setattr(namespace, self.dest, [fn_outer(i) for i in filenames])
        pass
    pass

class MergeBinaryBedFilesAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        sys.stderr.write("Merging %d bbg files\n" % len(filenames))
        bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = local_ref_version)
        for f in filenames:
            bt.read_as_regions(f)
            pass
        setattr(namespace, self.dest, bt)
        pass
    pass

class IntersectBinaryBedFilesAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        sys.stderr.write("intersecting %d bbg files\n" % len(filenames))
        bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = local_ref_version)
        bt.read_to_bases(myBedTools.binarybedfilegenome, filenames[0], myBedTools.set_to_one)
        for f in filenames[1:]:
            bt.read_to_bases(myBedTools.binarybedfilegenome, f, myBedTools.bitfn_and)
            gc.collect()
            pass
        setattr(namespace, self.dest, bt)
        pass
    pass


def munge_regions(opts):
    if opts.regions != None: sys.stderr.write('Initial number of bases considered (regions): %s\n' % locale.format('%d', opts.regions.total_length(), grouping=True))
    
    if opts.regions != None and opts.intersect_region != None:
        opts.regions.bases &= opts.intersect_region.bases
        opts.intersect_region = None
        pass

    if opts.regions != None and opts.exclude_region != None:
        sys.stderr.write('Excluding %s bases.\n' % locale.format('%d', opts.exclude_region.total_length(), grouping=True))
        opts.regions.bases &= ~opts.exclude_region.bases
        opts.exclude_region = None

    elif opts.regions == None and opts.exclude_region != None:
        sys.stderr.write('Excluding %s bases.\n' % locale.format('%d', opts.exclude_region.total_length(), grouping=True))
        opts.regions = opts.exclude_region
        opts.exclude_region = None
        opts.regions.bases.invert()
        pass
        
    if opts.regions != None: sys.stderr.write('Final number of bases considered (regions, exclude, intersect, etc): %s\n\n' % locale.format('%d', opts.regions.total_length(), grouping=True))

    if opts.regions == None and not opts.vcf_is_ms_file:
        sys.stderr.write('NO MASKING GIVEN - ASSUMING WHOLE GENOME IS PERFECT (hint: this is probably not a good assumption)\n\n')
        pass

    return


class ancestral_vcf(object):

    def __init__(self, filename):
        
        self.vcf = {}
        vcffile = open(filename, 'r')

        c = 0
        
        for line in vcffile:

            if line.strip().startswith('#'): continue

            try:
                [chrom, pos, _, ref, alt, qual, _, _, _, gt_info] = line.strip().split()
                pos = int(pos)
            except ValueError:
                print "Too many lines in ANCESTRAL VCF: %s?" % filename 
                print "Expecting one individual (10 columns):"
                print line
                sys.exit(-1)
                pass

            if chrom not in self.vcf:
                self.vcf[chrom] = self.base_chrom_dict()
                pass

            if pos in self.vcf[chrom]:
                print "error - duplicate position in VCF file?"
                print chrom, pos, ref, alt
                print line
                sys.exit(-1)
                pass

            self.add_site(chrom, pos, ref)
            
            c += 1
            pass

        return

    @staticmethod
    def base_chrom_dict():
        return {}

    def init_chrom(self, chrom):
        self.vcf[chrom] = vcf_class.base_chrom_dict()
        return

    def add_site(self, chrom, pos, ref):
        self.vcf[chrom][pos] = ref
        return

    def has_site(self, chrom, pos):
        return pos in self.vcf[chrom]

    def get_base_one_based(self, chrom, pos, illumina_chrs = False):
        if not self.has_site(chrom, pos):
            return 'N'
        return self.vcf[chrom][pos]
