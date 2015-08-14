import argparse, sys
from collections import defaultdict
import locale
import gzip


class vcf_class(object):

    geno_map = {'0/0':0, '0/1':1, '1/0':1, '1/1':2,
                '0|0':0, '0|1':1, '1|0':1, '1|1':2}

    def __init__(self, filename = None, ancestral_bsg = None, vcf_has_illumina_chrnums = False, regions = None):
        
        if filename == None:
            # blank file - just initialize the vcf_class
            vcffile = []
        elif filename.endswith('.gz'):
            vcffile = gzip.open(filename)
        else:
            vcffile = open(filename, 'r')
            pass
        
        self.vcf = {}
        if filename != None: sys.stderr.write("Reading VCF file %s..\n" % filename)
        c = 0
        warn_strip_chrom = False
        
        for line in vcffile:
            
            if line.strip().startswith('#'): continue

            try:
                [chrom, pos, _, ref, alt, qual, _, _, _, gt_info] = line.strip().split()
                pos = int(pos)
            except ValueError:
                print "Too many lines in ARCHAIC VCF: %s?" % filename 
                print "Expecting one individual (10 columns):"
                print line
                sys.exit(-1)
                pass

            #print opts
            #opts != None and 
            if vcf_has_illumina_chrnums and chrom.startswith('chr'):
                chrom = chrom[3:]
                warn_strip_chrom = True
                pass

            # opts != None and 
            if regions != None and not regions.in_region_one_based(chrom, pos):
                # if opts.debug: print "skipping archaic site [not in regions]", chrom, pos, filename
                continue

            if chrom not in self.vcf:
                # assume that a site that's not in the vcf does actually have the MH ref (we're filtering out non-callable sites, so that's ok)
                # missingdict allows us to not add something to the defaultdict just because we query it
                # IF THIS IS CHANGED, ALSO CHANGE IT IN READ_MS.PY! (should be set via a function..)
                self.init_chrom(chrom)
                pass

            if pos in self.vcf[chrom]:
                print "error - duplicate position in VCF file?"
                print chrom, pos, ref, alt
                print line
                sys.exit(-1)
                pass

            if len(alt) > 1:
                print "Error reading archaic VCF file:"
                print vcffile
                print "len(alt) > 1: require only bi-allelic SNPs in archaic VCF"
                print line
                sys.exit(-1)
                pass

            gt = gt_info[:3]
            # alts = [s.upper() for s in alt.split(',')]

            anc = None
            if ancestral_bsg != None:
                anc = ancestral_bsg.get_base_one_based(chrom, pos)
                pass

            self.add_site(chrom, pos, gt, ref, alt, anc)
            
            # print chrom, pos, ref, alt, qual, gt_info, self.vcf[chrom][pos]
            c += 1
            pass

        if filename != None: sys.stderr.write(" with %d lines.\n" % c)
        if warn_strip_chrom:
            sys.stderr.write(" WARNING: REMOVED LEADING 'chr' FROM CHROMOSOME NAMES, TO MATCH ILLUMINA VCFS\n")
            pass
        
        return

    @staticmethod
    def base_chrom_dict():
        # return missingdict(lambda : (True, '0/0', 'N', ['N']))
        return {}

    def init_chrom(self, chrom):
        self.vcf[chrom] = vcf_class.base_chrom_dict()
        return

    def add_site(self, chrom, pos, gt, ref, alt, ancestral):
        # old: IF THIS IS CHANGED, ALSO CHANGE IT IN READ_MS.PY! (should be set via a function..)
        # (has_ref, genotype ('0/0', etc), ref_base, alt_base) # NOT ALLOWED TO HAVE MULTIPLE ALTS

        if ancestral == alt:
            gt = gt.replace('1', '.').replace('0', '1').replace('.', '0')
            self.vcf[chrom][pos] = ('0' in gt, gt, ancestral, ref)
        else:
            self.vcf[chrom][pos] = ('0' in gt, gt, ref, alt)
            pass
        return

    def has_derived(self, chrom, pos):
        if not self.has_site(chrom, pos):
            return False
        return '1' in self.vcf[chrom][pos][1]

    def get_derived(self, chrom, pos):
        if not self.has_site(chrom, pos):
            return 'N'
        return self.vcf[chrom][pos][3]

    def get_derived_count(self, chrom, pos):
        if not self.has_site(chrom, pos):
            return 0
        return self.vcf[chrom][pos][1].count('1')

    def has_site(self, chrom, pos):
        return pos in self.vcf[chrom]

    def get_derived_sites(self, chrom, winstart, winend):
        return [p for p in xrange(winstart, winend) \
                    if self.has_derived(chrom, p)]

    def get_derived_sites_with_der_count(self, chrom, winstart, winend):
        return [(p, self.get_derived_count(chrom,p)) for p in xrange(winstart, winend) \
                    if self.has_derived(chrom, p)]
        

    # def has_var_at_site(self, chrom, pos):
    #     return pos in self.vcf[chrom]

    # def get_genotype_str_one_based(self, chrom, pos):
    #     return self.vcf[chrom][pos][1]

    # def get_genotype_one_based(self, chrom, pos, ref_is_derived = False):
    #     return (2-vcf_class.geno_map[self.vcf[chrom][pos][1]]) if ref_is_derived else vcf_class.geno_map[self.vcf[chrom][pos][1]]

    # def get_alts_one_based(self, chrom, pos):
    #     return self.vcf[chrom][pos][3]

    # def get_ref_one_based(self, chrom, pos):
    #     return self.vcf[chrom][pos][2]

    # def has_ref(self, chrom, pos):
    #     return self.vcf[chrom][pos][0]

    # def get_bases_one_based(self, chrom, pos):
    #     #print self.vcf[chrom][pos]
    #     #print self.get_alts_one_based(chrom, pos), '+', ([self.get_ref_one_based(chrom, pos)] if self.has_ref(chrom, pos) else [])
    #     return list(set(self.get_alts_one_based(chrom, pos) + ([self.get_ref_one_based(chrom, pos)] if self.has_ref(chrom, pos) else [])))

    # # add_chr_to_name is in case we are using vcf chromosomes (but the bsg requires chr1, etc - this should be fixed)
    # def has_derived(self, chrom, pos, ancbsg, snp = None, add_chr_to_name = False):
    #     anc = ancbsg.get_base_one_based(chrom, pos, add_chr_to_name)
    #     if anc == 'N':
    #         return False

    #     ## we're not making sure it *shares* a derived with the snp, just checking to see if this has a derived allele
    #     if snp == None:
    #         arc_bases = [b.upper() for b in self.get_bases_one_based(chrom, pos)]
    #         if anc not in arc_bases:
    #             return True
    #         elif len(arc_bases) > 1:
    #             return True
    #         return False

    #     ## in this case it must *share* a derived with the snp
    #     else:

    #         # first make sure that the reference base on this snp is ancestral (it's supposed to be set like that in read_vcf)
    #         if snp['ref'].upper() != anc:
    #             print "Error: 'reference' in snp isn't ancestral, but it should be set that way in read_vcf."
    #             print "Ancestral:", anc
    #             print "SNP", snp
    #             sys.exit(-1)
    #             pass
            
    #         return snp['alt'].upper() in [b.upper() for b in self.get_bases_one_based(chrom, pos)]
        
    #     print "Error: shouldn't be able to get this far in function has_derived"
    #     return
        
    # def get_archaic_derived_positions(self, chrom, winstart, winend, ancbsg, regions = None):

    #     return [p for p in xrange(winstart, winend) \
    #                 if self.has_derived(chrom, p, ancbsg) and \
    #                 (regions == None or regions.in_region_one_based(full_chrom, p))]
    
    pass


def process_archaic_vcfs(archaic_vcf_filenames, archaic_bsg, opts):
    return [vcf_class(f, archaic_bsg, opts.vcf_has_illumina_chrnums, opts.regions) for f in archaic_vcf_filenames]

# class VCFFileAction(argparse.Action):
#     def __call__(self, parser, namespace, filelist, option_string=None):

#         vcf_list = list()

#         for filename in filelist:

        
#         setattr(namespace, self.dest, vcf_list)
#         # setattr(namespace, self.dest, vcf)
#         pass
#     pass
