from __future__ import division
import sys, traceback, os, gzip
from array import array
from bitarray import bitarray
import bitarray as bitarray_m
sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/mapped_snps/')
sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/metrics/')
from region_stats import region_type_stats
from BaseLookup import BaseLookup
import fileinput
from operator import itemgetter

class myBedTools:

    bedfile = 'b'
    binarybedfile = 'bb'
    binarybedfilegenome = 'bbg'
    binaryseqfilegenome = 'bsg'
    binaryseqfilegenomechr = 'bsgc'
    varfile = 'v'
    vcffile = 'vcf'
    vcffile_ref = 'vcfref'
    vcffile_var = 'vcfvar'
    primatesnpfile = 's'
    seattleseqfile = 'ss'
    ancestralstatefile = 'a'
    output_coverage = 'cov_out'
    basic_bases = set(['A', 'C', 'G', 'T'])

    bedtypes = (bedfile, binarybedfile, binarybedfilegenome)

    binarybedfilegenome_code = '0000000000000000000000000000000000000000000000000000000000101010'
    binaryseqfilegenome_code = '0000000000000000000000000000000000000000000000000000000000101011'

    binaryseq_decode = {'A': bitarray('001'),
                        'C': bitarray('010'),
                        'G': bitarray('011'),
                        'N': bitarray('000'),
                        'T': bitarray('100')}

    def __init__(self, output_type, debug = False, output_file_name = None, invert = False, ref_version = 'b37', initialize = True):

        # store type of analysis
        self.output_type = output_type
        self.debug = False
        self.invert = invert
        initval = 0
        if self.invert: initval = 1
        self.initialize = initialize
        self.factor = 1
        if output_type == self.binaryseqfilegenome: self.factor = 3

        if ref_version == 'b37':
            self.genome_len = BaseLookup.chr_lens['genome']
            self.chr_offset = BaseLookup.chr_offset
            self.chr_lens = BaseLookup.chr_lens
            self.chrs = BaseLookup.chrs
        elif ref_version == 'b36':
            self.genome_len = BaseLookup.chr_lens_36['genome']
            self.chr_offset = BaseLookup.chr_offset_36
            self.chr_lens = BaseLookup.chr_lens_36
            self.chrs = BaseLookup.chrs
        elif ref_version == 'tair10':
            self.genome_len = BaseLookup.chr_lens_tair10['genome']
            self.chr_offset = BaseLookup.chr_offset_tair10
            self.chr_lens = BaseLookup.chr_lens_tair10
            self.chrs = BaseLookup.chrs_tair10
            pass
        sys.stderr.write('Creating bitarray of length %d\n' % (self.genome_len * self.factor))
        self.bases = bitarray(self.genome_len * self.factor)
        sys.stderr.write('Initializing bitarray.. ')
        #if self.initialize: self.bases[:] = bool(initval)
        if self.initialize: self.bases.setall(initval)
        sys.stderr.write('Finished\n')
        self.output_file_name = output_file_name

        if output_type == self.output_coverage:
            self.bases = array('B', self.bases)
            self.max_cov = 255
            pass

        pass


    @staticmethod
    def open_file(filename, print_header = False, discard_header = False):
        if filename == 'stdin' or filename == '-':
            file = sys.stdin
        elif filename.endswith('z'):
            file = gzip.open(filename, 'r')
        else:
            file = open(filename, 'r')
            pass
        if print_header: print file.readline().strip()
        elif discard_header: file.readline().strip()
        return file

    ## the return_base functionality for the various file types is really nonstandard!
    def parse_line(self, filetype, line, return_base = False):
        #if line.strip().startswith('#'): self.parse_line(filetype, line, return_base = return_base)
        if filetype == self.bedfile:
            [c,s,e] = line.strip().split()[:3]
            s = int(s)
            e = int(e)
            #print'%s %d' % (c,s)
            if return_base:
                a = line.strip().split()[3]
                # if it's anything weird ('AAA' for an indel), then set it to just missing
                if a not in self.basic_bases: a = 'N'
                return (c, s, e, a)
            pass
        elif filetype == self.varfile:
            [_,c,s,e] = line.strip().split()[:4]
            s = int(s)
            e = int(e)
        elif filetype == self.vcffile or filetype == self.vcffile_ref or filetype == self.vcffile_var:
            [c,e,_,r,a] = line.strip().split()[:5]
            e = int(e)
            s = e-1
            # if the alt allele is just ., then it is actually the ref allele
            if a == '.' and not filetype == self.vcffile_var: a = r
            # if it's anything weird ('AAA' for an indel), then set it to just missing
            if a not in self.basic_bases: a = 'N'
            # return alt
            if return_base and filetype == self.vcffile: return (c, s, e, a)
            if return_base and filetype == self.vcffile_var: return (c, s, e, a)
            if return_base and filetype == self.vcffile_ref: return (c, s, e, r)
            pass
        elif filetype == self.primatesnpfile:
            [c,e] = line.strip().split()[:2]
            e = int(e)
            s = e-1
            if return_base: return (c, s, e, line.strip().split()[3])
            pass
        elif filetype == self.ancestralstatefile:
            [_,c,e] = line.strip().split()[:3]
            e = int(e)
            s = e-1
            pass
        elif filetype == self.seattleseqfile:
            [_,c,e] = line.strip().split()[:3]
            c = 'chr' + c
            e = int(e)
            s = e-1
            pass
        else:
            print "error in parse line: bad filetype!"
            print filetype
            sys.exit(-1)
            pass
        return (c, s, e)

    # functions

    bitfn_and = 'bitfn_and'

    @staticmethod
    def set_to_one(state):
        return 1
    
    @staticmethod
    def set_to_zero(state):
        return 0
    
    @staticmethod
    def increment(state):
        return state+1

    @staticmethod
    def output_if_zero(state):
        return state == 0

    @staticmethod
    def output_if_one(state):
        return state == 1

    @staticmethod
    def output_if_two(state):
        return state == 2

    filename = 'none'
    def read_to_bases(self, filetype, filename, fn, exp_chr = None, header = False):
        self.filename = os.path.basename(filename)
        sys.stderr.write('Reading file...')
        if filetype == self.binarybedfile:
            if exp_chr == None:
                print "must send exp_chr for bb files!"
                sys.exit()
                pass
            ar = array('B')
            try:
                ar.fromfile(open(filename, 'rb'), self.genome_len+1)
            except EOFError:
                pass
            if len(ar) != self.chr_lens[exp_chr] and not (exp_chr == 'chrM' and len(ar) == self.chr_lens[exp_chr]+1):
                print "expected length of array to match given start and end!"
                print filename
                print exp_chr
                print len(ar), self.chr_lens[exp_chr]
                sys.exit(-1)
                pass
            for i in xrange(self.chr_lens[exp_chr]):
                if ar[i] == 1: self.bases[self.chr_offset[exp_chr]+i] = fn(self.bases[self.chr_offset[exp_chr]+i])
                pass
            pass
        
        elif filetype == self.binarybedfilegenome:
            ar = bitarray()
            code = bitarray()
            file = open(filename, 'r')
            #code.fromfile(file)
            code.fromfile(file, 8)
            if code.to01() != myBedTools.binarybedfilegenome_code:
                print "unexpected code for binary bed file genome!"
                print "code:     ", code
                print "expected: ", myBedTools.binarybedfilegenome_code
                sys.exit(-1)
            try:
                #ar.fromfile(file, self.genome_len + 16)
                ar.fromfile(file)
            except EOFError:
                pass
            if ar.length() != bitarray_m.bits2bytes(self.genome_len) * 8:
                print "expected length of array to match given start and end!"
                print filename
                print ar.length(), bitarray_m.bits2bytes(self.genome_len), self.genome_len
                sys.exit(-1)
                pass
            # shorten ar to genome_len
            for i in xrange(ar.length() - self.genome_len):
                #print self.genome_len, ar.length(), 'popping'
                ar.pop()
                pass
            #### COULD CHECK FNS FOR BITARRAY STUFF (AND, OR ETC)
            if not self.initialize:
                self.initialize = True
                if fn == myBedTools.set_to_one:
                    self.bases = ar
                elif fn == myBedTools.set_to_zero:
                    self.bases = ~ar
                else:
                    print "not initialized, and using fn other than set_to_one!"
                    print fn
                    print len(self.bases)
                    sys.exit(-1)
                    return
                pass
            elif fn == myBedTools.set_to_one:
                self.bases |= ar
            elif fn == myBedTools.set_to_zero:
                self.bases &= ~ar
            elif fn == myBedTools.bitfn_and:
                self.bases &= ar
            else:
                for i in xrange(self.end):
                    if ar[i] == 1: self.bases[i] = fn(self.bases[i])
                    pass
                pass
            pass

        elif filetype == self.binaryseqfilegenome:
            ar = bitarray()
            code = bitarray()
            if filename.endswith('z'):
                file = gzip.open(filename, 'r')
            else:
                file = open(filename, 'r')
                pass
            code.fromfile(file, 8)
            if code.to01() != myBedTools.binaryseqfilegenome_code:
                print "unexpected code for binary seq file genome!"
                print "code:     ", code
                print "expected: ", myBedTools.binaryseqfilegenome_code
                sys.exit(-1)
            try:
                #ar.fromfile(file, self.genome_len * self.factor + 16)
                ar.fromfile(file)
            except EOFError:
                pass
            if ar.length() != bitarray_m.bits2bytes(self.genome_len * self.factor) * 8:
                print "expected length of array to match given start and end!"
                print filename
                print ar.length(), bitarray_m.bits2bytes(self.genome_len * self.factor), self.genome_len * self.factor
                sys.exit(-1)
                pass
            # shorten ar to genome_len
            for i in xrange(ar.length() - self.genome_len * self.factor):
                #print self.genome_len, ar.length(), 'popping'
                ar.pop()
                pass

            self.bases = ar
            pass

        elif filetype == self.binaryseqfilegenomechr:
            print "THIS DOESN'T WORK"
            sys.exit(-1)
            ar = bitarray()
            code = bitarray()
            file = open(filename, 'r')
            code.fromfile(file, 8)
            if code.to01() != myBedTools.binaryseqfilegenome_code:
                print "unexpected code for binary seq file genome!"
                print "code:     ", code
                print "expected: ", myBedTools.binaryseqfilegenome_code
                sys.exit(-1)
            try:
                chr_byte_start = self.chr_offset[exp_chr] * self.factor // 8
                chr_partial_byte_start = (self.chr_offset[exp_chr] * self.factor) % 8
                print 'starting at', chr_byte_start, 'bytes, with', chr_partial_byte_start, 'offset'
                file.seek(chr_byte_start, 0)
                ar.fromfile(file, self.chr_lens[exp_chr] * self.factor + 16)
                print ar[0:100]
            except EOFError:
                pass
            if ar.length() != bitarray_m.bits2bytes(self.chr_lens[exp_chr] * self.factor) * 8:
                print "expected length of array to match given start and end!"
                print filename
                print ar.length(), bitarray_m.bits2bytes(self.genome_len * self.factor), self.genome_len * self.factor
                sys.exit(-1)
                pass
            # shorten ar to genome_len
            for i in xrange(ar.length() - self.genome_len * self.factor):
                #print self.genome_len, ar.length(), 'popping'
                ar.pop()
                pass

            self.bases = ar
            pass

        elif self.output_type == myBedTools.binaryseqfilegenome:
            infile = myBedTools.open_file(filename, discard_header = header)
            
            for line in infile:
                if line.strip().startswith('#'): continue
                #print line
                [l_chr, l_start, l_end, l_base] = self.parse_line(filetype, line, return_base = True)

                if len(l_base) != 1: continue
                
                if self.debug: print 'bed', l_chr, l_start, l_end, l_base
                site = self.chr_offset[l_chr] * self.factor + l_start * self.factor
                self.bases[site : site + self.factor] = self.binaryseq_decode[l_base.upper()]
                pass
            pass
        else:
            infile = myBedTools.open_file(filename, discard_header = header)
            
            for line in infile:
                [l_chr, l_start, l_end] = self.parse_line(filetype, line)
                
                if self.debug: print 'bed', l_chr, l_start, l_end
                for i in xrange(l_start, l_end):
                    self.bases[self.chr_offset[l_chr]+i] = fn(self.bases[self.chr_offset[l_chr]+i])
                    pass
                pass
            pass
        sys.stderr.write(' done\n')
        return

    
    def read_and_output_file(self, filetype, filename, fn, header = False):
        sys.stderr.write('Reading file, outputing file\n')
        infile = myBedTools.open_file(filename, print_header = header)

        for line in infile:
            [l_chr, l_start, l_end] = self.parse_line(filetype, line)

            if self.debug: print 'bed', l_chr, l_start, l_end
            #### COULD CHANGE THIS TO PERCENT AND HAVE A "MAX TO PASS" FLAG, LIKE BEDOPS
            if self.amount_in_region(l_chr, l_start, l_end) > 0 and fn(1):
                print line.strip('\n')
                pass
            pass
        return




    #### REALLY NEED TO FIX THIS!  DOES NOT WORK AT ALL RIGHT NOW
    #### SWITCH TO NP.ARRAY(BITARRAY(CHRLEN), NP.INT8) FOR EACH CHR.  MAKE SURE THAT WE CAN EASILY SAVE NP.ARRAY TO FILE!
    #### NP.ARRAY HAS THE ADVANTAGE OF BEING ABLE TO DO SELF.BASES[CHR] += AR[CHR_OFFSET[CHR]:CHR_OFFSET[CHR]+CHR_LENS[CHR]], TO GET COVERAGE (UP TO 127, WITH NO WARNING)
    #### can save to file using np.save('filename', data) and np.load('filename')
    def report_coverage_for_windows(self, filetype, filename, outputtype = 'hist', individuals = None):
        infile = myBedTools.open_file(filename)

        for line in infile:
            [l_chr, l_start, l_end] = self.parse_line(filetype, line)

            if self.debug: print 'window', l_start, l_end, max(l_start, self.start), min(l_end, self.end)

            if l_chr != self.chr or l_start < self.start or l_end > self.end:
                print "window outside declared range!"
                print 'window', l_start, l_end, max(l_start, self.start), min(l_end, self.end)
                sys.exit(-1)
                pass

            if outputtype == 'hist':
                h = []
                for i in xrange(l_start, l_end):
                    c = self.bases[i-self.start]
                    if self.debug: print ' ', i, i-self.start, c, len(h)
                    if c >= len(h):
                        x = [0] * (c+1)
                        x[:len(h)] = h
                        h = x
                        pass
                    h[c] += 1
                    pass
                for i in xrange(len(h)):
                    if h[i] > 0: print "%s\t%d\t%d\t%d\t%d" % (l_chr, l_start, l_end, i, h[i])
                    pass
                pass
            if outputtype == 'histdropone':
                h = [[0] * (len(individuals)+1) for i in xrange(len(individuals)+1)]
                for i in xrange(l_start, l_end):
                    ind_cov = [self.get_drop_one_coverage(ind, self.bases[i-self.start]) for ind in individuals]
                    tot_cov = sum(ind_cov)
                    drop_one_cov = [tot_cov - c for c in ind_cov]
                    if self.debug: print ' ', i, i-self.start, tot_cov, ind_cov, drop_one_cov, h
                    h[tot_cov][0] += 1
                    for ind in individuals:
                        h[drop_one_cov[ind]][ind+1] += 1
                        pass
                    pass
                for i in xrange(len(h)):
                    output = "%s\t%d\t%d\t%d\t%d" % (l_chr, l_start, l_end, i, h[i][0])
                    for ind in individuals:
                        output += "\t%d" % h[i][ind+1]
                        pass
                    print output
                    pass
                pass
            pass
        return


    def bases_to_bed(self, fn):

        sys.stderr.write('Bases to bed\n')
        
        if self.output_type == self.binarybedfilegenome:
            file = open(self.output_file_name, 'wb')
            code = bitarray(myBedTools.binarybedfilegenome_code)
            code.tofile(file)
            self.bases.tofile(file)
            file.close()
            
        elif self.output_type == self.binaryseqfilegenome and self.output_file_name == None:
            for c in self.chrs:
                # for c in ['chr%s' % c for c in [13, 14, 15, 20, 21, 22, 'X']]:
                s = self.chr_offset[c]
                e = s + self.chr_lens[c]
                sys.stderr.write('outputing bed for %s, %d-%d\n' % (c, s, e))
                for i in xrange(s,e):
                #for i in xrange(s+85232,s+1085332):
                    # print c,i,s,e
                    b = self.get_base_one_based(c, i-s+1)
                    if b in self.basic_bases:
                        print '%s\t%d\t%d\t%s' % (c, i-s, i-s+1, b)
                        pass
                    pass
                pass
            
        elif self.output_type == self.binaryseqfilegenome:
            file = open(self.output_file_name, 'wb')
            code = bitarray(myBedTools.binaryseqfilegenome_code)
            code.tofile(file)
            self.bases.tofile(file)
            file.close()
            
        else:

            ohone = bitarray('01')
            for c in self.chrs:
                s = self.chr_offset[c]
                e = s + self.chr_lens[c]
                sys.stderr.write('outputing bed for %s, %d-%d\n' % (c, s, e))
                chr_bases = self.bases[s:e]
                pre = []
                if self.bases[s] == 1:
                    pre.append(0)
                    pass
                for i in pre + chr_bases.search(ohone):
                    j = i+2
                    while chr_bases[j]:
                        j += 1
                        continue
                    print '%s\t%d\t%d' % (c, i+1, j)
                    pass
                pass
            pass
        
        return

    def pad_beds_by(self, n):

        ohone = bitarray('01')
        oneoh = bitarray('10')

        for c in self.chrs:
            s = self.chr_offset[c]
            e = s + self.chr_lens[c]
            sys.stderr.write('padding bed for %s, %d-%d\n' % (c, s, e))
            chr_bases = self.bases[s:e]
            chr_len = self.chr_lens[c]
            pre = []
            if chr_bases[0] == 1: pre.append(0)
            start_positions = pre + chr_bases.search(ohone)
            post = []
            if chr_bases[-1] == 1: post.append(len(chr_bases)-1)
            end_positions = chr_bases.search(oneoh) + post
            if len(start_positions) != len(end_positions):
                print "error in start / end positions"
                print len(start_positions), len(end_positions)
                pass
            for index, bed_start in enumerate(start_positions):
                bed_end = end_positions[index]
                j = bed_end+1
                #if bed_start > 2244000: print c, bed_start+1, j, max(0,bed_start+1-n), bed_start+1, j, min(j+n,chr_len)
                chr_bases[max(0,bed_start+1-n) : bed_start+1] = True
                chr_bases[j : min(j+n,chr_len)] = True
                #print '%s\t%d\t%d' % (c, i+1, j)
                pass
            self.bases[s:e] = chr_bases
            pass
        return

    
    ### region fns
    ### ONLY ALLOWS ONE REGIONS FILE TO BE READ!!!

    def check_pos(self, chr, start, end, fail = True):
        if start < 0 or end > self.chr_lens[chr]:
            if fail:
                print "in myBedTools, position outside genome def!"
                print 'pos', chr, start, end, self.chr_lens[chr]
                traceback.print_exc(file=sys.stdout)
                sys.exit(-1)
            else:
                sys.stderr.write( "in myBedTools, position outside genome def!\n" )
                sys.stderr.write( ' '.join(str(s) for s in ('pos', chr, start, end, self.chr_lens[chr], '\n') ) )
                return False
            pass
        return True

    def read_as_regions(self, filename, file_type = binarybedfilegenome, header = False, invert = False):
        
        if invert or self.invert:
            fn = self.set_to_zero
        else:
            fn = self.set_to_one
        
        self.read_to_bases(file_type, filename, fn, header = header)
        return
    
    def get_base_one_based(self, chr, pos):
        self.check_pos(chr, pos-1, pos)
        site = self.chr_offset[chr] * self.factor + pos * self.factor
        return self.bases[site - self.factor : site].decode(self.binaryseq_decode)[0]
    
    def get_base_zero_based(self, chr, pos):
        self.check_pos(chr, pos, pos+1)
        site = self.chr_offset[chr] * self.factor + pos * self.factor
        return self.bases[site : site + self.factor].decode(self.binaryseq_decode)[0]

    def total_length(self):
        return self.bases.count()

    def in_region_one_based(self, chr, pos):
        self.check_pos(chr, pos-1, pos)
        return self.bases[self.chr_offset[chr] + pos - 1] > 0
    
    def in_region_zero_based(self, chr, pos):
        self.check_pos(chr, pos, pos+1)
        return self.bases[self.chr_offset[chr] + pos] > 0
    
    def amount_in_region(self, chr, start, end, fail = True):
        
        ## if the full region isn't part of the chr definition, then either
        ##   a) if fail==True, check_pos will sys.exit()
        ##   b) if fail==False, check_pos will return False, and we'll adjust the size of the region and return the number of mapped bases

        if self.check_pos(chr, start, end, fail):
            return sum(self.bases[self.chr_offset[chr]+start : self.chr_offset[chr]+end])
        
        if start > self.chr_lens[chr]: return 0
        if end > self.chr_lens[chr]: end = self.chr_lens[chr]
        
        return sum(self.bases[self.chr_offset[chr]+start : self.chr_offset[chr]+end])
    
    def percent_in_region(self, chr, start, end):
        return self.amount_in_region(chr, start, end) / (end-start)
        
    pass






if __name__ == "__main__":

    
    if len(sys.argv) < 2:
        print 'useage:'
        print '%s analysistype [-b bedfile] [-bb binarybedfile chr] [-bbg binarybedfile] [-bsg binaryseqfile] [-v varfile] [-vcf vcffile] [-vcfref vcffile (use ref allele)] [-s primatesnpfile] [-a ancestralstatefile] [-ss seattleseqfile] [-obbg filename] [-obsg filename]' % sys.argv[0]
        print "  bedfile can be '-' or 'stdin' to use stdin."
        sys.exit()
        pass

    analysis_type = sys.argv[1]

    output_type = myBedTools.bedfile
    output_file_name = None
    ref_version = 'b37'

    # get file data
    files = []
    i = 2
    filetypes = {'-b': myBedTools.bedfile,
                 '-bb': myBedTools.binarybedfile,
                 '-bbg': myBedTools.binarybedfilegenome,
                 '-bsg': myBedTools.binaryseqfilegenome,
                 '-bsgc': myBedTools.binaryseqfilegenomechr,
                 '-v': myBedTools.varfile,
                 '-vcf': myBedTools.vcffile,
                 '-vcfref': myBedTools.vcffile_ref,
                 '-vcfvar': myBedTools.vcffile_var,
                 '-s': myBedTools.primatesnpfile,
                 '-ss': myBedTools.seattleseqfile,
                 '-a': myBedTools.ancestralstatefile}
    while i < len(sys.argv):
        flag = sys.argv[i]
        if flag in filetypes:
            ifn = sys.argv[i+1]
            f = {}
            f['name'] = ifn
            f['type'] = filetypes[flag]
            files.append(f)
            i += 2
            if f['type'] == myBedTools.binarybedfile or f['type'] == myBedTools.binaryseqfilegenomechr:
                f['chr'] = sys.argv[i]
                i += 1
                pass
            elif flag == '-bsg':
                output_type = myBedTools.binaryseqfilegenome
                output_file_name = None
                pass
            pass
        elif flag[0] != '-':
            ifn = sys.argv[i]
            f = {}
            f['name'] = ifn
            f['type'] = myBedTools.binarybedfilegenome
            files.append(f)
            i += 1
        elif flag == '-h':
            files[-1]['header'] = True
            i += 1
            pass
        elif flag == '-obbg':
            output_type = myBedTools.binarybedfilegenome
            output_file_name = sys.argv[i+1]
            i += 2
            pass
        elif flag == '-obsg':
            output_type = myBedTools.binaryseqfilegenome
            output_file_name = sys.argv[i+1]
            i += 2
            pass
        elif flag == '-b36':
            ref_version = 'b36'
            i += 1
            pass
        elif flag == '-tair10':
            ref_version = 'tair10'
            i += 1
            pass
        pass

    # should we initialize?
    initialize = True
    # first file is binary bed
    # first file is binary seq
    # filtering based on binary bed
    if files[0]['type'] == myBedTools.binarybedfilegenome \
           or files[0]['type'] == myBedTools.binaryseqfilegenome \
           or files[0]['type'] == myBedTools.binaryseqfilegenomechr \
           or (analysis_type == 'filter' and files[1]['type'] == myBedTools.binarybedfilegenome):
        initialize = False
        pass

    # make bedtools object
    bto = myBedTools(output_type, output_file_name = output_file_name, ref_version = ref_version, initialize = initialize)

    # perform analysis!
    if analysis_type == 'report':
        f = files[0]
        bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_one, header = 'header' in f)
        bto.bases_to_bed(myBedTools.output_if_one)

    elif analysis_type == 'merge':
        for f in files:
            exp_chr = None
            if f['type'] == myBedTools.binarybedfile: exp_chr = f['chr']
            if f['type'] == myBedTools.binaryseqfilegenomechr: exp_chr = f['chr']
            bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_one, exp_chr = exp_chr, header = 'header' in f)
            pass
        bto.bases_to_bed(myBedTools.output_if_one)

    elif analysis_type == 'invert':
        f = files[0]
        bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_zero, header = 'header' in f)
        bto.bases_to_bed(myBedTools.output_if_one)

    elif analysis_type == 'create-bsg':
        if files[0]['type'] != myBedTools.primatesnpfile:
            print "bsg creation requires one primatesnpfile"
            print files[0]['type']
            sys.exit(-1)
            pass
        f = files[0]
        bto.read_to_bases(f['type'], f['name'], None)
        bto.bases_to_bed(None)

    elif analysis_type == 'subtract':

        print >> sys.stderr, "subtract seems to be broken!!"
        print "subtract seems to be broken!!"
        sys.exit(-1)
        
        if files[0]['type'] in myBedTools.bedtypes:
            f = files[0]
            bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_one)
            for f in files[1:]:
                bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_zero)
                pass
            bto.bases_to_bed(myBedTools.output_if_one)
        else:
            for f in files[1:]:
                bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_one)
                pass
            f = files[0]
            bto.read_and_output_file(f['type'], f['name'], myBedTools.output_if_zero, header = 'header' in f)
            pass
        pass
    
    elif analysis_type == 'intersect':
        f = files[0]
        bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_one)
        for f in files[1:]:
            bto.read_to_bases(f['type'], f['name'], myBedTools.bitfn_and)
            pass
        bto.bases_to_bed(myBedTools.output_if_one)

    elif analysis_type == 'length':
        for f in files:
            bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_one)
            pass
        print bto.total_length()

    elif analysis_type == 'filter':
        f = files[1]
        bto.read_to_bases(f['type'], f['name'], myBedTools.set_to_one, header = 'header' in f)
        f = files[0]
        bto.read_and_output_file(f['type'], f['name'], myBedTools.output_if_one, header = 'header' in f)
        pass

    elif analysis_type == 'coverage':
        window_file = files.pop()
        for f in files:
            bto.read_to_bases(f['type'], f['name'], myBedTools.increment)
            pass
        bto.report_coverage_for_windows(window_file['type'], window_file['name'])
        pass

    elif analysis_type == 'coverage-dropone':
        window_file = files.pop()
        for i,f in enumerate(files):
            bto.read_to_bases(f['type'], f['name'], myBedTools.drop_one_coverage(i))
            pass
        bto.report_coverage_for_windows(window_file['type'], window_file['name'], outputtype = 'histdropone', individuals = range(len(files)))
        pass

    elif analysis_type == 'regiontest':
        bto.read_as_regions('/net/akey/vol1/home/bvernot/tishkoff/primate_sequences/alignments/hg19.panTro2.synNet.axt.5.unmapped.zerobased.filtered.cpg.chr11', invert = True)
        print bto.amount_in_region('chr11', 100000, 110000)
        for i in range(100001, 110000):
            print i, bto.in_region('chr11', i)
            pass
        pass

    pass
    



