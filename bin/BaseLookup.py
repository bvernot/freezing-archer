from __future__ import division
import sys, os



######
##  one-based lookup

class BaseLookup:

    def __init__(self, basedir = './', ref_version = 'b37'):
        self.basedir = basedir
        if not basedir.endswith('/'):
            self.basedir += '/'
            pass
        self.seq_files = {}
        self.seq_files_header_len = {}
        self.bases_per_line = {}

        self.chr_len_hash = {}
        self.chr_len_hash['hg19'] = self.chr_lens
        self.chr_len_hash['hg18'] = self.chr_lens_36
        self.chr_len_hash['tair10'] = self.chr_lens_tair10
        
        if ref_version == 'b36':
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
        
        pass

    
    def getBase(self, species, chr, site):

        if site < 1 or site > self.chr_len_hash[species][chr]:
            print "bad base pos in BaseLookup:", chr, site
            print "too low or too high"
            sys.exit(-1)
            pass

        if not species in self.seq_files:
            self.seq_files[species] = {}
            self.seq_files_header_len[species] = {}
            self.bases_per_line[species] = {}
            pass

        if not chr in self.seq_files[species]:
            d = self.basedir + species + '/' + chr + '.fa'
            self.seq_files[species][chr] = open(d, 'r')
            self.seq_files_header_len[species][chr] = len(self.seq_files[species][chr].readline())
            self.bases_per_line[species][chr] = len(self.seq_files[species][chr].readline().strip())
            pass

        site_trans = self.seq_files_header_len[species][chr]
        site_trans += site - 1
        site_trans += (site-1) // self.bases_per_line[species][chr]
        
        self.seq_files[species][chr].seek(site_trans)
        b = self.seq_files[species][chr].read(1)
        if b == '\n':
            print "bad base pos (base returned is newline):", chr, site
            print 'next pos', self.seq_files[species][chr].read(1)
            print 'next pos', self.seq_files[species][chr].read(1)
            print 'next pos', self.seq_files[species][chr].read(1)
            self.seq_files[species][chr].seek(site_trans-1)
            print 'prev pos', self.seq_files[species][chr].read(1)
            self.seq_files[species][chr].seek(site_trans-2)
            print 'prev pos-1', self.seq_files[species][chr].read(1)
            sys.exit(-1)
        return b

    chr_nums = {'sim' : 1,
                'chr1' : 1,
                'chr2' : 2,
                'chr3' : 3,
                'chr4' : 4,
                'chr5' : 5,
                'chr6' : 6,
                'chr7' : 7,
                'chr8' : 8,
                'chr9' : 9,
                'chr10' : 10,
                'chr11' : 11,
                'chr12' : 12,
                'chr13' : 13,
                'chr14' : 14,
                'chr15' : 15,
                'chr16' : 16,
                'chr17' : 17,
                'chr18' : 18,
                'chr19' : 19,
                'chr20' : 20,
                'chr21' : 21,
                'chr22' : 22,
                'chrX'  : 23,
                'chrY'  : 24,
                'chrM'  : 25}

    chrs = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX' , 'chrY' , 'chrM')

    chrX_par1 = 2699520 ## end of par1
    chrX_par2 = 154931044 ## beginning of par2
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

    chr_lens_36 = {'chr1' : 247249719,
                   'chr2' : 242951149,
                   'chr3' : 199501827,
                   'chr4' : 191273063,
                   'chr5' : 180857866,
                   'chr6' : 170899992,
                   'chr7' : 158821424,
                   'chr8' : 146274826,
                   'chr9' : 140273252,
                   'chr10' : 135374737,
                   'chr11' : 134452384,
                   'chr12' : 132349534,
                   'chr13' : 114142980,
                   'chr14' : 106368585,
                   'chr15' : 100338915,
                   'chr16' : 88827254,
                   'chr17' : 78774742,
                   'chr18' : 76117153,
                   'chr19' : 63811651,
                   'chr20' : 62435964,
                   'chr21' : 46944323,
                   'chr22' : 49691432,
                   'chrX' : 154913754,
                   'chrY' : 57772954,
                   'chrM' : 16571}
    #'chrM' : 16570}
    
    
    chr_lens['genome'] = sum([chr_lens[c] for c in chrs])
    chr_lens_36['genome'] = sum([chr_lens_36[c] for c in chrs])

    chr_offset = {}
    running_sum = 0
    chr_offset_36 = {}
    running_sum_36 = 0
    for i,c in enumerate(chrs):
        chr_offset[c] = running_sum
        running_sum += chr_lens[c]
        chr_offset_36[c] = running_sum_36
        running_sum_36 += chr_lens_36[c]
        pass


    chr_lens_tair10 = {'chr1' : 30427671,
                   'chr2' : 19698289,
                   'chr3' : 23459830,
                   'chr4' : 18585056,
                   'chr5' : 26975502,
                   'chrC' : 154478,
                   'chrM' : 366924}
    chrs_tair10 = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chrC' , 'chrM')
    chr_nums_tair10 = {c:i for i,c in enumerate(chrs_tair10)}
    chr_lens_tair10['genome'] = sum([chr_lens_tair10[c] for c in chrs_tair10])

    chr_offset_tair10 = {}
    running_sum_tair10 = 0
    for i,c in enumerate(chrs_tair10):
        chr_offset_tair10[c] = running_sum_tair10
        running_sum_tair10 += chr_lens_tair10[c]
        pass


#import BaseLookup
if __name__ == "__main__":
    b = BaseLookup('/net/akey/vol1/home/bvernot/tishkoff/primate_sequences/')
    for i in range(1,51):
        print b.getBase('hg19', 'chrTest', i)
        pass

    for i in range(51,60):
        print b.getBase('hg19', 'chrTest', i)
        pass

    print '----'
    for i in range(100,155):
        print b.getBase('hg19', 'chrTest', i)
        pass
    
