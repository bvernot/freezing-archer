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

    
    pop_mapping = {'GS19700-1100-37-ASM' : 'asw1',
                   'GS19701-1100-37-ASM' : 'asw2',
                   'GS19703-1100-37-ASM' : 'asw3',
                   'GS19704-1100-37-ASM' : 'asw4',
                   'GS19834-1100-37-ASM' : 'asw5',
                   'GS06985-1100-37-ASM' : 'ceu1',
                   'GS06994-1100-37-ASM' : 'ceu2',
                   'GS07357-1100-37-ASM' : 'ceu3',
                   'GS10851-1100-37-ASM' : 'ceu4',
                   'GS12004-1100-37-ASM' : 'ceu5',
                   'GS12877-1100-37-ASM' : 'ceu6',
                   'GS12878-1100-37-ASM' : 'ceu7',
                   'GS12879-1100-37-ASM' : 'ceu8',
                   'GS12880-1100-37-ASM' : 'ceu9',
                   'GS12881-1100-37-ASM' : 'ceu10',
                   'GS12882-1100-37-ASM' : 'ceu11',
                   'GS12883-1100-37-ASM' : 'ceu12',
                   'GS12884-1100-37-ASM' : 'ceu12',
                   'GS12885-1100-37-ASM' : 'ceu13',
                   'GS12886-1100-37-ASM' : 'ceu14',
                   'GS12887-1100-37-ASM' : 'ceu15',
                   'GS12888-1100-37-ASM' : 'ceu16',
                   'GS12889-1100-37-ASM' : 'ceu17',
                   'GS12890-1100-37-ASM' : 'ceu18',
                   'GS12891-1100-37-ASM' : 'ceu19',
                   'GS12892-1100-37-ASM' : 'ceu20',
                   'GS12893-1100-37-ASM' : 'ceu21',
                   'GS18526-1100-37-ASM' : 'chb1',
                   'GS18537-1100-37-ASM' : 'chb2',
                   'GS18555-1100-37-ASM' : 'chb3',
                   'GS18558-1100-37-ASM' : 'chb4',
                   'GS20845-1100-37-ASM' : 'gih1',
                   'GS20846-1100-37-ASM' : 'gih2',
                   'GS20847-1100-37-ASM' : 'gih3',
                   'GS20850-1100-37-ASM' : 'gih4',
                   'GS000001807-ASM' : 'hdz1',
                   'GS000003124-ASM' : 'hdz2',
                   'GS000003125-ASM' : 'hdz3',
                   'GS00319-DNA_A02_1100_37-ASM' : 'hdz4',
                   'GS00319-DNA_B02_1100_37-ASM' : 'hdz5',
                   'GS18940-1100-37-ASM' : 'jpt1',
                   'GS18942-1100-37-ASM' : 'jpt2',
                   'GS18947-1100-37-ASM' : 'jpt3',
                   'GS18956-1100-37-ASM' : 'jpt4',
                   'GS19017-1100-37-ASM' : 'lwk1',
                   'GS19020-1100-37-ASM' : 'lwk2',
                   'GS19025-1100-37-ASM' : 'lwk3',
                   'GS19026-1100-37-ASM' : 'lwk4',
                   'GS21732-1100-37-ASM' : 'mkk1',
                   'GS21733-1100-37-ASM' : 'mkk2',
                   'GS21737-1100-37-ASM' : 'mkk3',
                   'GS21767-1100-37-ASM' : 'mkk4',
                   'GS19648-1100-37-ASM' : 'mxl1',
                   'GS19649-1100-37-ASM' : 'mxl2',
                   'GS19669-1100-37-ASM' : 'mxl3',
                   'GS19670-1100-37-ASM' : 'mxl4',
                   'GS19735-1100-37-ASM' : 'mxl5',
                   'HG00731-1100-37-ASM' : 'pur1',
                   'HG00732-1100-37-ASM' : 'pur2',
                   'HG00733-1100-37-ASM' : 'pur3',
                   'GS00319-DNA_A01_1100_37-ASM' : 'pyg1',
                   'GS00319-DNA_B01_1100_37-ASM' : 'pyg2',
                   'GS00319-DNA_C01_1100_37-ASM' : 'pyg3',
                   'GS00319-DNA_D01_1100_37-ASM' : 'pyg4',
                   'GS00319-DNA_E01_1100_37-ASM' : 'pyg5',
                   'GS000001577-ASM' : 'san1',
                   'GS00319-DNA_C02_1100_37-ASM' : 'san2',
                   'GS00319-DNA_D02_1100_37-ASM' : 'san3',
                   'GS00319-DNA_E02_1100_37-ASM' : 'san4',
                   'GS00319-DNA_G02_1100_37-ASM' : 'san5',
                   'GS20502-1100-37-ASM' : 'tsi1',
                   'GS20509-1100-37-ASM' : 'tsi2',
                   'GS20510-1100-37-ASM' : 'tsi3',
                   'GS20511-1100-37-ASM' : 'tsi4',
                   'GS18501-1100-37-ASM' : 'yri1',
                   'GS18502-1100-37-ASM' : 'yri2',
                   'GS18504-1100-37-ASM' : 'yri3',
                   'GS18505-1100-37-ASM' : 'yri4',
                   'GS18508-1100-37-ASM' : 'yri5',
                   'GS18517-1100-37-ASM' : 'yri6',
                   'GS19129-1100-37-ASM' : 'yri7',
                   'GS19238-1100-37-ASM' : 'yri8',
                   'GS19239-1100-37-ASM' : 'yri9'}

    pop_mapping_general = dict([(i, pop_mapping[i][:3]) for i in pop_mapping])

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
    
