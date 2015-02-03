import sys, itertools, operator


# ms 42 1 -s 10 
# 44126 40565 42561

# //
# segsites: 10
# positions: 0.1717 0.2230 0.2277 0.4523 0.4598 0.5201 0.7094 0.8533 0.9100 0.9894 
# 0010000001
# 0110000001
# 1000000000
# 1000000000
# 0010000001


#@profile
def process_ms_block_to_genotypes(ms_file , opts):


    ##  assume that the // line has just been consumed

    ## get segsites
    line = ms_file.readline().strip()
    while not line.startswith('segsites'):
        line = ms_file.readline().strip()
        if line == '':
            print "EOF"
            return None
        pass

    segsites = int(line.split()[1])

    ## get position lists
    while not line.startswith('positions'):
        line = ms_file.readline().strip()
        pass
    
    positions = list(itertools.imap(float, line.split()[1:]))
    print positions



    ## read genotypes, one ind at a time
    in_genotypes = False
    ind_genotypes = [None] * opts.num_inds
    for ind_num in range(opts.num_inds):

        c1_line = ms_file.readline().strip()
        ## make sure we're actually at a genotype line (could be an L or ages line)
        while not in_genotypes and c1_line[0] not in ('0', '1'):
            print "SKIPPING LINE:", c1_line
            c1_line = ms_file.readline().strip()
            pass
        in_genotypes = True
        c2_line = ms_file.readline().strip()
        print c1_line
        print c2_line

        ## check to make sure that we actually have two genotype lines
        if c1_line[0] not in ('0', '1') or c2_line[0] not in ('0', '1'):
            print "bad genotype line?"
            print c1_line
            print c2_line
            sys.exit(-1)
            pass
        
        c1_line = itertools.imap(int, c1_line)
        c2_line = itertools.imap(int, c2_line)

        ind_genotypes[ind_num] = itertools.imap(operator.add, c1_line, c2_line)
        pass
    
    snps = ['hey'] * segsites

    for snp_num, genotypes in enumerate(itertools.izip(*ind_genotypes)):
        snp_d = {}
        snp_d['pos'] = int(positions[snp_num] * opts.window_length)

        t_gt = list(itertools.imap(lambda j : genotypes[j], opts.updated_ind_indexed_to_orig_file_target))
        r_gt = list(itertools.imap(lambda j : genotypes[j], opts.updated_ind_indexed_to_orig_file_reference))
        print t_gt
        print r_gt
        genotypes = t_gt + r_gt

        snp_d['genotypes'] = genotypes
        snp_d['target'] = sum(t_gt) > 0
        snp_d['reference'] = sum(r_gt) > 0

        snps[snp_num] = snp_d
        
        print snp_d
        print opts
        pass



    print snps

    return snps
    
    # snp_d = {}
    
    # split_line = line.strip().split()
    #     #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096 HG00097
    # snp_d['chrom'] = split_line[0]
    # snp_d['pos'] = int(split_line[1])
    # snp_d['ref'] = split_line[3]
    # snp_d['alt'] = split_line[4]
    # snp_d['qual'] = split_line[5]
    # snp_d['filt'] = split_line[6]
    # snp_d['info'] = split_line[7]
    # format = split_line[8]

    # ## remove non-biallelic sites
    # if len(snp_d['ref']) != 1 or len(snp_d['alt']) != 1:
    #     if opts.debug: print "dropping snp because it is not biallelic, or is an indel"
    #     return None

    # ## get genotypes
    # genotypes = split_line[9:]

    # ## remove individual genotypes that aren't in the target or reference, and reorder (target first, then ref)
    # ## (reordered ind ids in read_1kg_ind_pop_file)
    # r_gt = list(itertools.imap(lambda j : genotypes[j], opts.updated_ind_indexed_to_orig_file_reference))
    # t_gt = list(itertools.imap(lambda j : genotypes[j], opts.updated_ind_indexed_to_orig_file_target))


    # genotypes = t_gt + r_gt
    # #print t_gt
    # #print r_gt
    # #print genotypes
    # t_gt = set(t_gt)
    # r_gt = set(r_gt)

    # snp_d['target'] = '0|1' in t_gt or '1|0' in t_gt or '1|1' in t_gt
    # snp_d['reference'] = not('0|0' in r_gt and len(r_gt) == 1)

    # ## finally, remove variants that aren't in the target or reference
    # if not snp_d['target'] and not snp_d['reference']:
    #     if opts.debug: print "dropping snp because it's NOT present in either the TARGET or REFERENCE"
    #     return None


    # ## get the location of genotypes, if necessary
    # if True or ':' in format:
    #     gt_loc = format.split(':')
    #     if 'GT' not in gt_loc:
    #         print "bad format?  GT not found"
    #         print format
    #         print split_line
    #         sys.exit(-1)
    #         pass
    #     gt_loc = gt_loc.index('GT')
    
    #     ## strip out just genotypes
    #     # genotypes[:] = [gt.split(':')[gt_loc] for gt in genotypes]
    #     genotypes[:] = list(itertools.imap(lambda gt : gt.split(':')[gt_loc], genotypes))
    #     pass
    
    # if opts.debug: print 'reading', snp_d['pos'], genotypes

    # gt_map = {'0|0':0, '1|0':1, '0|1':1, '1|1':2}
    # gt2 = list(itertools.imap(lambda gt : gt_map[gt], genotypes))
    # genotypes[:] = list(itertools.imap(lambda gt : gt_map[gt], genotypes))
    # if opts.debug: print genotypes


    # snp_d['genotypes'] = genotypes



    # if opts.debug: print

    # return snp_d



#@profile
def ms_to_genotypes_windowed(ms_file, winlen, winstep, ms_ind_pop_file, opts, start = 0):

    ## get sample ID order from the header (ignore comments) / results are saved to opts
    read_ms_header(ms_file, opts)

    ## process sample ID info / results are saved to opts
    read_ms_ind_pop_file(ms_ind_pop_file, opts)


    snps = []
    keep_reading = True
    chrom = 1
    winstart = start
    winend = winlen

    while keep_reading:

        snps = process_ms_block_to_genotypes(ms_file, opts)
        
        if snps == None: break
        
        yield (chrom, winstart, winend, snps)
        chrom += 1
        pass

    pass





def read_ms_ind_pop_file(f, opts):

    # sample  pop     super_pop       gender
    # HG00096 GBR     EUR     male
    # HG00097 GBR     EUR     female
    # HG00099 GBR     EUR     female

    header = f.readline()
    opts.sample_to_pop = {}
    opts.sample_to_superpop = {}
    opts.all_pops = set()
    opts.all_superpops = set()
    opts.updated_reference_population_inds = set()
    opts.updated_target_population_inds = set()

    for line in f:
        (sample, pop, superpop, _) = line.strip().split()
        opts.sample_to_pop[sample] = pop
        opts.sample_to_superpop[sample] = superpop
        opts.all_pops.add(pop)
        opts.all_superpops.add(superpop)

        if pop      in opts.reference_populations: opts.updated_reference_population_inds.add(sample)
        if superpop in opts.reference_populations: opts.updated_reference_population_inds.add(sample)
        if pop      in opts.target_populations:    opts.updated_target_population_inds.add(sample)
        if superpop in opts.target_populations:    opts.updated_target_population_inds.add(sample)

        pass

    opts.updated_target_population_inds.update(opts.target_individuals)
    opts.updated_reference_population_inds.update(opts.reference_individuals)

    opts.updated_ind_ids_target    = [i for i in opts.original_file_ind_ids if i in opts.updated_target_population_inds]
    opts.updated_ind_ids_reference = [i for i in opts.original_file_ind_ids if i in opts.updated_reference_population_inds]
    opts.updated_ind_indexed_to_orig_file_target    = [opts.original_file_ind_ids.index(ind) for ind in opts.updated_ind_ids_target]
    opts.updated_ind_indexed_to_orig_file_reference = [opts.original_file_ind_ids.index(ind) for ind in opts.updated_ind_ids_reference]

    opts.updated_ind_ids = opts.updated_ind_ids_target + opts.updated_ind_ids_reference
    # opts.updated_ind_indexed_to_orig_file = opts.updated_ind_indexed_to_orig_file_target + opts.updated_ind_indexed_to_orig_file_reference

    
    return# (sample_to_pop, sample_to_superpop)


# ms 42 1 -s 10 
# 44126 40565 42561

# //
# segsites: 10
# positions: 0.1717 0.2230 0.2277 0.4523 0.4598 0.5201 0.7094 0.8533 0.9100 0.9894 
# 0010000001
# 0110000001
# 1000000000
# 1000000000
# 0010000001

def read_ms_header(ms_file, opts):

    opts.ms_command = ms_file.readline().strip().split()
    opts.ms_seeds = ms_file.readline().strip().split()

    opts.num_inds = int(opts.ms_command[1]) / 2

    line = ms_file.readline()
    while line != '' and not line.startswith('//'):
        line = ms_file.readline()
        pass

    opts.original_file_ind_ids = ['i%d' % i for i in range(opts.num_inds)]
    #opts.original_file_ind_index = {id:i for i,id in enumerate(opts.original_file_ind_ids)}
    return
