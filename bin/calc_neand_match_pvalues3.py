from __future__ import division
import sys, os, argparse, numpy
from time import clock, time
from collections import defaultdict, deque
sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/filter_files/')
from myBedTools3 import myBedTools
from operator import itemgetter
from math import log, ceil
from itertools import islice
import locale
import random
import subprocess
locale.setlocale(locale.LC_ALL, 'en_US')
import tempfile
import gzip
import tables

CURRENT_FORMAT_VERSION = 4 # changed nhaps from uint8 to int16

start_time = time()
current_time = time()


parser = argparse.ArgumentParser(description='Manage haplotype data for entire genome, for the purpose of calculating pvalues for putative neanderthal haplotypes.')

# build parsers for various arguments
base_parser = argparse.ArgumentParser(add_help = False)
base_parser.add_argument('-d', '--debug', required=False, action='store_true')
base_parser.add_argument('-dp', '--debug-pval', required=False, action='store_true')
base_parser.add_argument('-chunk', '--chunk-to-process', required=False, default=None, type=int, help="the portion of the haplotype file to compute pvalues for.  STARTS AT 1")
base_parser.add_argument('-chunksize', '--chunksize', required=False, default=1000, type=int)
base_parser.add_argument('-o', '--output-file', required=False, type=argparse.FileType('w'), default=None)
base_parser.add_argument('-prog', '--print-progress', default = 0, const = 100000, type = int, nargs = '?')
base_parser.add_argument('-l', '--limit', default = None, type = int)

h5_parser = argparse.ArgumentParser(add_help = False)
h5_parser.add_argument('-df', '--database-file', required=True)

h5_dest_parser = argparse.ArgumentParser(add_help = False)
h5_dest_parser.add_argument('-odf', '--output-database-file', required=True)

hapfiles_parser = argparse.ArgumentParser(add_help = False)
hapfiles_parser.add_argument('-f', '--haplotype-table-file', required=True)

tolerance_parser = argparse.ArgumentParser(add_help = False)
tolerance_parser.add_argument('-len-tol', '--len-tolerance', required=False, default=1000, type=int)
tolerance_parser.add_argument('-freq-tol', '--freq-tolerance', required=False, default=.1, type=float)
tolerance_parser.add_argument('-sd-tol', '--sd-tolerance', required=False, default=.1, type=float)
tolerance_parser.add_argument('-sites-tol', '--sites-tolerance', required=False, default=1, type=int)

putative_regions_parser = argparse.ArgumentParser(add_help = False)
putative_regions_parser.add_argument('-pred', '--predicted-regions-file', required=True, type=argparse.FileType('r'))

compression_parser = argparse.ArgumentParser(add_help = False)
compression_parser.add_argument('-c', '--compression-algorithm', required=False, default=None)
compression_parser.add_argument('-clev', '--compression-level', required=False, type = int, default=1)

expected_rows_parser = argparse.ArgumentParser(add_help = False)
expected_rows_parser.add_argument('-er', '--expected-rows', required=False, type=int, default = 500000000)

remove_regions_parser = argparse.ArgumentParser(add_help = False)
remove_regions_parser.add_argument('-rm', '--remove-regions', required=False, nargs=3, action='append', default=[])

subparsers = parser.add_subparsers(dest='func')
create_database_parser = subparsers.add_parser('create', parents = [base_parser, h5_parser, hapfiles_parser, compression_parser, expected_rows_parser, remove_regions_parser])
test_database_parser = subparsers.add_parser('test', parents = [base_parser, h5_parser])
copy_database_parser = subparsers.add_parser('copy', parents = [base_parser, h5_parser, h5_dest_parser, compression_parser, expected_rows_parser])
get_pvals_parser = subparsers.add_parser('pvals', parents = [base_parser, h5_parser, putative_regions_parser, tolerance_parser])

opts = parser.parse_args()

if opts.output_file != None:
    sys.stdout = opts.output_file
    pass

if opts.debug: print opts
if opts.debug: print

if opts.chunk_to_process != None and opts.chunk_to_process < 1:
    print "chunk to process cannot be less than one!"
    print opts.chunk_to_process
    sys.exit(-1)
    pass

class Haplotype(tables.IsDescription):

    # names = ['count', 'mapped_bases_bin', 'len', 'mh_sites', 'tot_sites', 'sfs', 'std_dev', 'match'])
    # 1 10000 11111 20 28 141 62 14.0
    # 1 10000 11111 20 32 106 82 10.0
    # 6 10000 11111 21 36 170 51 15.0

    ## OLD
    # nhaps, chr, start, stop, len, mh_sites, mh_alt, mh_n_ref, all_sites, all_alt, freq, std_dev
    # 2       1       1002932 1007222 4290    9       5       5       10      5       0.522280856054  0.253164851538
    # 1       1       1014864 1019180 4316    9       6       3       11      6       0.4957490472    0.250605448462

    count = tables.Int16Col(pos=0)
    mapped_bases_bin = tables.UInt16Col(pos=1)
    hap_len = tables.UInt16Col(pos=2)
    mh_sites = tables.Int16Col(pos=3)
    tot_sites = tables.Int16Col(pos=4)
    sfs = tables.Int16Col(pos=5)
    std_dev = tables.Int16Col(pos=6)
    match = tables.Float32Col(pos=7)
    pass


if opts.func == 'create':

    opts.haplotype_table_file = gzip.open(opts.haplotype_table_file)

    mydata = tables.openFile(opts.database_file, 'w', 'haplotypes')
    
    if opts.compression_algorithm == None:
        filters = None
    else:
        filters = tables.Filters(complevel=opts.compression_level, complib=opts.compression_algorithm)
        pass

    hap_table = mydata.createTable(mydata.root, "table", Haplotype, filters = filters, expectedrows = opts.expected_rows)
    hap_rows = hap_table.row
    hap_table.setAttr('format_version', CURRENT_FORMAT_VERSION)

    print "inserting rows"

    for linenum, line in enumerate(opts.haplotype_table_file):
        
        if opts.print_progress > 0 and linenum % opts.print_progress == 0:
            sys.stderr.write('inserting line %d: %s\n' % (linenum, line))
            pass

        if opts.limit != None and linenum > opts.limit:
            sys.stderr.write('reached limit at %d: %s\n' % (linenum, line))
            break
        
        line = line.strip().split()

        if len(line) != 8: 
            print "bad line", line, linenum
            continue
        
        # names = ['count', 'mapped_bases_bin', 'len', 'mh_sites', 'tot_sites', 'sfs', 'std_dev', 'match'])
        # 1 10000 11111 20 28 141 62 14.0
        # 1 10000 11111 20 32 106 82 10.0
        # 6 10000 11111 21 36 170 51 15.0
        
        count, mapped_bases_bin, hap_len, mh_sites, tot_sites, sfs, std_dev, match = line
        
        # filter to: len > 10000 & mapped_bases_bin > 5000 & tot_sites >= 8?
        
        # if sum([line[1] == c and int(line[2]) <= int(e) and int(s) <= int(line[3]) for (c,s,e) in opts.remove_regions]) > 0:
        #     # print "line overlaps banned region", line, linenum, [line[1] == c and int(line[2]) <= int(e) and int(s) <= int(line[3]) for (c,s,e) in opts.remove_regions]
        #     continue
            
        # nhaps, chr, start, stop, len, mh_sites, mh_alt, mh_n_ref, all_sites, all_alt, freq, std_dev
        # 2       1       1002932 1007222 4290    9       5       5       10      5       0.522280856054  0.253164851538
        # 1       1       1014864 1019180 4316    9       6       3       11      6       0.4957490472    0.250605448462

        #try:
            
        hap_rows['count'] = count
        hap_rows['mapped_bases_bin'] = mapped_bases_bin
        hap_rows['hap_len'] = hap_len
        hap_rows['mh_sites'] = mh_sites
        hap_rows['tot_sites'] = tot_sites
        hap_rows['sfs'] = sfs
        hap_rows['std_dev'] = std_dev
        hap_rows['match'] = match
        hap_rows.append()
            
        # except Exception as e:
        #     print "error when inserting row"
        #     print line
        #     print linenum
        #     raise e
        pass
    
    print 'processed rows', opts.haplotype_table_file, time() - current_time

    current_time = time()

    hap_table.flush()

    # current_time = time()
    # hap_table.cols.mh_sites.createIndex(0, 'ultralight')
    # hap_table.flush()
    # print 'indexed', time()-current_time

    current_time = time()
    hap_table.cols.hap_len.createIndex(0, 'ultralight')
    hap_table.flush()
    print 'indexed', time()-current_time

    # current_time = time()
    # hap_table.cols.mapped_bases_bin.createIndex(0, 'ultralight')
    # hap_table.flush()
    # print 'indexed', time()-current_time

    current_time = time()
    a = [ x[:] for x in hap_table.where('(hap_len > 25000) & (hap_len < 30000) & (mh_sites > 7) & (mh_sites < 10) & (std_dev > 4) & (std_dev < 15)')]
    print mh_sites, len(a), time() - current_time

    pass

    # print 'tables', hap_tables

    # mydata.root.haplotypes.cols.hap_len.createCSIndex()
    # mydata.root.haplotypes.cols.freq.createCSIndex()
    # mydata.root.haplotypes.cols.std_dev.createCSIndex()
    # mydata.root.haplotypes.cols.mh_sites.createCSIndex()

    # print 'finished indexing', time() - current_time
    # current_time = time()

    mydata.close()
    sys.exit()
    pass

elif opts.func == 'copy':

    mydata = tables.openFile(opts.database_file, 'r', 'haplotypes')
    mydata_new = tables.openFile(opts.output_database_file, 'w', 'haplotypes')
    
    if opts.compression_algorithm == None:
        filters = None
    else:
        filters = tables.Filters(complevel=opts.compression_level, complib=opts.compression_algorithm)
        pass
    hap_table = mydata.createTable(mydata.root, "haplotypes", Haplotype, filters = filters, expectedrows = opts.expected_rows)
    hap_table.setAttr('format_version', CURRENT_FORMAT_VERSION)
    hap_row = hap_table.row


elif opts.func == 'test':

    mydata = tables.openFile(opts.database_file, 'r', 'haplotypes')
    #print mydata
    #print repr(mydata)

    hap_table = mydata.root.table

    if hap_table.getAttr('format_version') != CURRENT_FORMAT_VERSION:
        print "incorrect format version for this dataset:"
        print "file version:", hap_table.getAttr('format_version')
        print "expected version:", CURRENT_FORMAT_VERSION
        sys.exit(-1)
        pass

    # current_time = time()
    # a = [ x[:] for x in hap_table.where('(hap_len > 25000) & (hap_len < 30000) & (mapped_bases_bin > 15000) & (mapped_bases_bin < 20000)')]
    # print 'read', time()-current_time, len(a)

    random.seed(12345)

    for _ in range(2):
        current_time = time()
        a = [ x[:] for x in hap_table.where('(hap_len > 25000) & (hap_len < 30000) & (mapped_bases_bin > 15000) & (mapped_bases_bin < 20000) & ' +
                                            '(mh_sites > 20) & (mh_sites < 30)')]
        print 'read a', time()-current_time, len(a)

        current_time = time()
        a = [ x[:] for x in hap_table.where('(hap_len >= 25000) & (hap_len <= 27000) & (mapped_bases_bin >= 15000) & (mapped_bases_bin <= 17000) & ' +
                                            '(mh_sites >= 25) & (mh_sites <= 27)')]
        print 'read b', time()-current_time, len(a)
        pass

    start_rand_time = time()
    for _ in range(10):
        current_time = time()
        l=random.randint(20000,45000)
        m=random.randint(int(l*.5),int(l*.75))
        mh=random.randint(20,70)
        a = [ x[:] for x in hap_table.where('(hap_len > %d - 1000) & (hap_len < %d + 1000) & ' % (l,l) + 
                                            '(mapped_bases_bin >= %d - 1000) & (mapped_bases_bin <= %d + 1000) & ' % (m, m) +
                                            '(mh_sites >= %d-1) & (mh_sites <= %d+1)' % (mh,mh))]
        print 'read r:', l, m, mh, time()-current_time, len(a)
        pass
    print 'TOTAL RAND TIME:', time()-start_rand_time

    # print "IT DOES NOT WORK TO COMPARE LIKE THIS: 30000 > hap_len > 25000"
    # current_time = time()
    # a = [ x[:] for x in hap_table.where('(30000 > hap_len > 25000) & (mapped_bases_bin > 15000) & (mapped_bases_bin < 20000) & ' +
    #                                     '(mh_sites > 20) & (mh_sites < 30)')]
    # print 'read', time()-current_time, len(a)


#    print "pval"

    # names = ['count', 'mapped_bases_bin', 'len', 'mh_sites', 'tot_sites', 'sfs', 'std_dev', 'match'])
    # 1 10000 11111 20 28 141 62 14.0
    # 1 10000 11111 20 32 106 82 10.0
    # 6 10000 11111 21 36 170 51 15.0
    
    # num_haplotypes = sum(x[0] for x in a)
    # match_pct = .5
    # for match_pct in (x/10 for x in range(10)):
    #     match_pct_by_rows = [x[0] if x[7]/x[4] >= match_pct else 0 for x in a]
    #     pval = (sum(match_pct_by_rows) + 1) / (num_haplotypes + 1)
    #     print match_pct, pval, num_haplotypes
    #     pass


    # current_time = time()
    # a = [ x[:] for x in hap_table.where('(hap_len > 25000) & (hap_len < 30000) & (mapped_bases_bin > 15000) & (mapped_bases_bin < 20000)')]
    # print 'read', time()-current_time, len(a)

    # a = [ x[:] for x in hap_table.where('(hap_len < 30000) & (freq > .38) & (freq < .4) & (std_dev > .3) & (std_dev < .4)')]
    # a = [ x[:] for x in hap_table.where('(hap_len > 7000) & (hap_len < 10000) & (freq > .6) & (freq < .8) & (std_dev > .1) & (std_dev < .2)')]
    #print len(a)
    #print a
    #print hap_table[:10]

    pass




first_line = True
def report_line(line, first_line, random_haplotype = ['NA'] * 12, num_similar_haplotypes = 0, num_independent_haplotypes = 0, pvals = ['NA']*6, tol_factor = 1):
    if first_line:
        print '\t'.join(expected_header + ['similar_haplotypes', 'independent_haplotypes', "pval_alt", "pval_ref", "pval_alt_all", "pval_ref_all", "pval_rand_alt", "pval_rand_alt_all", 'tol_factor'] + ['rand_chr', 'rand_start', 'rand_stop', 'rand_mh_sites', 'rand_mh_alt', 'rand_mh_n_ref', 'rand_all_sites', 'rand_all_alt', 'rand_freq', 'rand_stdfreq'])
        first_line = False
        pass
    
    # list(putative_haplotype[5:]) + 
    print '\t'.join([str(s) for s in line + [num_similar_haplotypes, num_independent_haplotypes] + pvals + [tol_factor] + [random_haplotype[i] for i in (1,2,3,5,6,7,8,9,10,11)]])
    #print '\t'.join(line + random_haplotype + [str(s) for s in [num_similar_haplotypes, num_independent_haplotypes] + pvals + [tol_factor]])
    if opts.debug: print
    sys.stdout.flush()
    return first_line


if opts.func == 'pvals':

    mydata = tables.openFile(opts.database_file, 'r', 'haplotypes')
    
    #if opts.debug: print "haplotype table:", mydata.root.haplotypes, time() - current_time
    current_time = time()

    # import pandas
    # def load_hdf(s):
    #     t = time()
    #     if opts.debug: print "loading db into pandas", s
    #     a = pandas.read_hdf(opts.database_file, 'haplotypes_%d' % s)
    #     if opts.debug: print ".. db loading done:", time()-t
    #     return a

    # hap_tables = {load_hdf}
    

    header = opts.predicted_regions_file.readline().strip()

    expected_header = ['chrom', 'hap_start', 'hap_end', 'ind', 'start', 'stop', 'pval', 'num_tag_snps', 'num_snps', 'recomb', 'run_tag', 'subset_tag', 'tag_snps', 'neand_snps', 
                       'correct_hap', 'pct_tag_snps_on_hap', 'tag_snps_haps', 'num_neand_callable_tag_snps', 'callable_pos', 'callable_neand', 'callable_haps', 
                       'num_haplotype_sites', 'real_hap_start', 'real_hap_stop', 'num_runs', 'num_runs_bidir', 'num_match_in_subset', 'num_match_in_subset_bidir',
                       'tag_freq', 'tag_stdfreq',
                       'hap_len', 'mh_sites', 'mh_alt', 'mh_n_ref', 'all_sites', 'all_alt', 'freq', 'stdfreq']

    if header.split() != expected_header:
        print "header doesn't match expected header!"
        print header.split()
        print expected_header
        sys.exit(-1)
        pass
    
    num_haplotype_sites_col = expected_header.index('num_haplotype_sites')
    real_hap_start_col = expected_header.index('real_hap_start')
    real_hap_stop_col = expected_header.index('real_hap_stop')
    ind_col = expected_header.index('ind')
    correct_hap_col = expected_header.index('correct_hap')
    run_tag_col = expected_header.index('run_tag')
    # _col = expected_header.index('')


    for line_num, line in enumerate(opts.predicted_regions_file):

        if opts.chunk_to_process != None and line_num < (opts.chunk_to_process-1) * opts.chunksize: continue
        if opts.chunk_to_process != None and line_num >= opts.chunk_to_process * opts.chunksize: sys.exit(0)

        if line.startswith('file'): continue

        line = line.strip().split()

        chrom = line[0]
        chrom_num = chrom[3:]
        num_haplotype_sites = int(line[num_haplotype_sites_col])

        if num_haplotype_sites < 2:
            first_line = report_line(line, first_line)
            continue
        
        real_hap_start = int(line[real_hap_start_col])
        real_hap_stop = int(line[real_hap_stop_col])
        correct_hap = int(line[correct_hap_col])
        run_tag = line[run_tag_col]
        ind = int(line[ind_col])
         # = int(line[_col])
         # = int(line[_col])

        if correct_hap not in (1,2):
            if opts.debug: print "do not have correct haplotype:", correct_hap
            first_line = report_line(line, first_line)
            continue

        # nhaps, chr, start, stop, len, mh_sites, mh_alt, mh_n_ref, all_sites, all_alt, freq, std_dev
        # 2       1       1002932 1007222 4290    9       5       5       10      5       0.522280856054  0.253164851538
        putative_haplotype = ['x', chrom_num, real_hap_start, real_hap_stop] + line[-8:]
        mh_sites = putative_haplotype[5]
        # print putative_haplotype
        if opts.debug: print "putative/test haplotype", putative_haplotype

        if mh_sites == 'NA':
            first_line = report_line(line, first_line)
            continue

        if opts.debug: print "loading appropriate haplotype table"
        # hap_table = mydata.getNode('/haplotypes_%s' % mh_sites) # tagged on mh_sites
        hap_tables = [mydata.getNode('/haplotypes_%s' % (int(mh_sites) + tol)) for tol in range(-opts.sites_tolerance, opts.sites_tolerance+1) if (int(mh_sites) + tol) > 8] # tagged on mh_sites

        if len(hap_tables) == 0:
            print 'no hap tables?', mh_sites
            first_line = report_line(line, first_line)
            continue

        #print mh_sites
        #print hap_tables

        if opts.debug: print 'format version', hap_tables[0].getAttr('format_version')
        if opts.debug: print

        if hap_tables[0].getAttr('format_version') != CURRENT_FORMAT_VERSION:
            print "incorrect format version for this dataset:"
            print "file version:", hap_tables[0].getAttr('format_version')
            print "expected version:", CURRENT_FORMAT_VERSION
            sys.exit(-1)
            pass


        independent_haplotypes = []
        random_hap = None
        prev_tol_factor = 0
        tol_factor = 1/1.5 # this is the value needed to *start* at 1 (after multiplying by 1.5
        line_current_time = time()
        while len(independent_haplotypes) < 40 and tol_factor < 2.5:

            tol_factor = tol_factor * 1.5

            # nhaps, chr, start, stop, len, mh_sites, mh_alt, mh_n_ref, all_sites, all_alt, freq, std_dev
            # 2       1       1002932 1007222 4290    9       5       5       10      5       0.522280856054  0.253164851538
            pytable_request = ('(hap_len > %d - %d) & (hap_len < %d + %d) & (freq > %f - %f) & (freq < %f + %f) & (std_dev > %f - %f) & (std_dev < %f + %f)' % 
                               (int(putative_haplotype[4]), opts.len_tolerance * tol_factor, int(putative_haplotype[4]), opts.len_tolerance * tol_factor, 
                                float(putative_haplotype[10]), opts.freq_tolerance * tol_factor, float(putative_haplotype[10]), opts.freq_tolerance * tol_factor, 
                                float(putative_haplotype[11]), opts.sd_tolerance * tol_factor, float(putative_haplotype[11]), opts.sd_tolerance * tol_factor))

            if opts.debug: print "pytable request (tol factor %f): %s" % (tol_factor, pytable_request)
            
            current_time = time()
            # all_similar_haplotypes = [ x[:] for x in hap_table.where(pytable_request) ]
            all_similar_haplotypes = []
            for hap_table in hap_tables:
                all_similar_haplotypes += [ x[:] for x in hap_table.where(pytable_request) ]
                pass
            num_similar_haplotypes = len(all_similar_haplotypes)

            if opts.debug: print "\nidentified %d haplotypes [currently have %d independent haps]" % (num_similar_haplotypes, len(independent_haplotypes)), time() - line_current_time, time() - current_time
            
            if len(all_similar_haplotypes) == 0: continue
            
            all_similar_haplotypes.sort(key=lambda x : (str(x[1]), x[2], x[3]))
            f = tempfile.NamedTemporaryFile(delete=True)

            # nhaps, chr, start, stop, len, mh_sites, mh_alt, mh_n_ref, all_sites, all_alt, freq, std_dev
            # 2       1       1002932 1007222 4290    9       5       5       10      5       0.522280856054  0.253164851538
            f.writelines(('\t'.join((str(hap[i]) for i in (1,2,3,0,4,5,6,7,8,9,10,11))) + '\n' for hap in all_similar_haplotypes))
            f.flush()
            if opts.debug: print "\nwrote similar haps to tempfile %s" % f.name, time() - current_time

            map_ids = subprocess.check_output("bedops --ec -m %s | awk 'BEGIN {OFS=\"\t\"} {print $0, NR}' | bedmap --ec --echo-map-id %s -" % (f.name, f.name), shell=True).split()
            if opts.debug: print "got map ids", len(map_ids), time() - current_time

            current_time = time()
            map_id_dict = defaultdict(lambda:[])
            for i,hap in enumerate(all_similar_haplotypes):
                # map_id_dict[map_ids[i]].append(hap)
                ## make sure to count each haplotype in a region the number of times that it was present in the original data (hap[0])
                map_id_dict[map_ids[i]].extend((hap,) * hap[0])
                # print '  ', map_ids[i], len(map_id_dict[map_ids[i]]), hap
                pass
            if opts.debug: print "created full list of similar haps by map id", len(map_ids), time() - current_time
            # print "seeded map_id_dict", time() - current_time
            # print random.choice(map_id_dict['400']), time() - current_time
            # print random.choice(map_id_dict['401']), time() - current_time
            # print random.choice(map_id_dict['402']), time() - current_time

            #for map_id in map_id_dict:
            #    print map_id, len(map_id_dict[map_id]), map_id_dict[map_id]
            #    pass

            independent_haplotypes = [random.choice(map_id_dict[map_id]) for map_id in map_id_dict]
            if opts.debug: print "generated ind haps", len(independent_haplotypes), time() - current_time

            if random_hap == None:
                ## CHANGED BACK TO THIS METHOD OF SELECTING A RANDOM HAP, TO SEE IF IT AFFECTED THE RANDOM HAP PVAL DISTRIBUTION
                rand_id = random.choice(map_id_dict.keys())
                random_hap = random.choice(map_id_dict[rand_id])
                #random_hap = random.choice([item for sublist in map_id_dict.values() for item in sublist])
                if opts.debug: print "selected similar random haplotype:", random_hap
                pass

            pass
        
        if len(independent_haplotypes) == 0:
            if opts.debug: print "NO SIMILAR HAPLOTYPES FOUND"
            first_line = report_line(line, first_line)
            continue

        num_independent_haplotypes = len(independent_haplotypes)
        if opts.debug: print "identified %d independent haplotypes" % num_independent_haplotypes, time() - line_current_time
        current_time = time()
        
        
        # match_pct_alt = int(putative_haplotype[7]) / int(putative_haplotype[6])
        # match_pct_ref = int(putative_haplotype[8]) / int(putative_haplotype[6])
        # match_pct_alt_all = int(putative_haplotype[10]) / int(putative_haplotype[9])
        # match_pct_ref_all = int(putative_haplotype[8]) / int(putative_haplotype[9])
        
        ###### instead of doing each pct/pval independently, just streamline the code (they're in the same order as shown above)
        # nhaps, chr, start, stop, len, mh_sites, mh_alt, mh_n_ref, all_sites, all_alt, freq, std_dev
        # 2       1       1002932 1007222 4290    9       5       5       10      5       0.522280856054  0.253164851538
        pval_indices = ((6,5), (7,5), (9,8), (7,8), (6,5), (9,8))
        test_haplotypes = (putative_haplotype, putative_haplotype, putative_haplotype, putative_haplotype, random_hap, random_hap)
        
        pvals = []
        for i, (num,den) in enumerate(pval_indices):
            match_pct = int(test_haplotypes[i][num]) / int(test_haplotypes[i][den])
            match_pct_by_rows = [x[num]/x[den] >= match_pct for x in independent_haplotypes]
            pval = (sum(match_pct_by_rows) + 1) / (num_independent_haplotypes + 1)
            pvals.append(pval)
            if opts.debug: print "pval (%d,%d): %f" % (num, den, pval)
            pass
        
        first_line = report_line(line, first_line, random_hap, num_similar_haplotypes, num_independent_haplotypes, pvals, tol_factor)
        
        pass
    pass
