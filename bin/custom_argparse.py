import argparse, sys
from collections import defaultdict
import locale
import gzip

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

class vcf_class(object):

    geno_map = {'0/0':0, '0/1':1, '1/0':1, '1/1':2, 
                '0|0':0, '0|1':1, '1|0':1, '1|1':2}

    def __init__(self, vcf):
        self.vcf = vcf
        return

    def has_var_at_site(self, chrom, pos):
        return pos in self.vcf[chrom]

    def get_genotype_str_one_based(self, chrom, pos):
        return self.vcf[chrom][pos][1]

    def get_genotype_one_based(self, chrom, pos, ref_is_derived = False):
        return (2-vcf_class.geno_map[self.vcf[chrom][pos][1]]) if ref_is_derived else vcf_class.geno_map[self.vcf[chrom][pos][1]]

    def get_alts_one_based(self, chrom, pos):
        return self.vcf[chrom][pos][3]

    def get_ref_one_based(self, chrom, pos):
        return self.vcf[chrom][pos][2]

    def has_ref(self, chrom, pos):
        return self.vcf[chrom][pos][0]

    def get_bases_one_based(self, chrom, pos):
        #print self.vcf[chrom][pos]
        #print self.get_alts_one_based(chrom, pos), '+', ([self.get_ref_one_based(chrom, pos)] if self.has_ref(chrom, pos) else [])
        return self.get_alts_one_based(chrom, pos) + ([self.get_ref_one_based(chrom, pos)] if self.has_ref(chrom, pos) else [])

    pass

class VCFFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        if filename.endswith('.gz'):
            vcffile = gzip.open(filename)
        else:
            vcffile = open(filename, 'r')
            pass

        #vcf = defaultdict(lambda : ['N'])
        vcf = {}
        sys.stderr.write("Reading VCF file %s..\n" % filename)
        c = 0
        warn_strip_chrom = False
        for line in vcffile:
            if line.strip().startswith('#'): continue
            try:
                [chrom, pos, _, ref, alt, qual, _, _, _, gt_info] = line.strip().split()
            except ValueError:
                print "Too many lines in VCF?  Expecting one individual (9 columns):"
                print line
                sys.exit(-1)
                pass
            if namespace.vcf_has_illumina_chrnums and chrom.startswith('chr'):
                chrom = chrom[3:]
                warn_strip_chrom = True
                pass
            if chrom not in vcf:
                # assume that a site that's not in the vcf does actually have the MH ref (we're filtering out non-callable sites, so that's ok)
                # missingdict allows us to not add something to the defaultdict just because we query it
                vcf[chrom] = missingdict(lambda : (True, '0/0', 'N', ['N']))
                pass
            if int(pos) in vcf[chrom]:
                print "error - duplicate position in VCF file?"
                print chrom, pos, ref, alt
                print line
                sys.exit(-1)
                pass
            gt = gt_info[:3]
            # (has_ref, genotype ('0/0', etc), ref_base, alt_bases)
            vcf[chrom][int(pos)] = ('0' in gt, gt, ref, alt.split(','))
            # print chrom, pos, ref, alt, qual, gt_info, vcf[chrom][int(pos)]
            c += 1
            pass
        sys.stderr.write(" with %d lines.\n" % c)
        if warn_strip_chrom:
            sys.stderr.write(" WARNING: REMOVED LEADING 'chr' FROM CHROMOSOME NAMES, TO MATCH ILLUMINA VCFS\n")
            pass
        
        setattr(namespace, self.dest, vcf_class(vcf))
        # setattr(namespace, self.dest, vcf)
        pass
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

    if opts.regions == None: 
        sys.stderr.write('NO MASKING GIVEN - ASSUMING WHOLE GENOME IS PERFECT (hint: this is probably not a good assumption)\n\n')
        pass

    return
