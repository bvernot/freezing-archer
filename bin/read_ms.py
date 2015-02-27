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
def process_ms_block_to_genotypes(ms_file , num_inds, window_length):
    """
    Process one ms block into positions and genotypes.
    Returns ((pos, gts), ..); a generator that returns consecutive positions and genotypes.
    """

    ##  assume that the // line has just been consumed

    ## get segsites
    line = ms_file.readline().strip()
    while not line.startswith('segsites'):
        line = ms_file.readline().strip()
        if line == '':
            return None
        pass

    segsites = int(line.split()[1])

    ## get position lists
    while not line.startswith('positions'):
        line = ms_file.readline().strip()
        pass
    
    # map the space [0,window_length)
    positions = list(itertools.imap(lambda s : int(float(s) * (window_length-1)), line.split()[1:]))

    ## read genotypes, one ind at a time
    in_genotypes = False
    ind_genotypes = [None] * num_inds
    for ind_num in range(num_inds):

        c1_line = ms_file.readline().strip()
        ## make sure we're actually at a genotype line (could be an L or ages line)
        while not in_genotypes and c1_line[0] not in ('0', '1'):
            print "SKIPPING LINE:", c1_line
            c1_line = ms_file.readline().strip()
            pass
        in_genotypes = True
        c2_line = ms_file.readline().strip()
        # print c1_line
        # print c2_line

        ## check to make sure that we actually have two genotype lines
        if c1_line[0] not in ('0', '1') or c2_line[0] not in ('0', '1'):
            print "bad genotype line?"
            print c1_line
            print c2_line
            sys.exit(-1)
            pass
        
        c1_line = itertools.imap(int, c1_line)
        c2_line = itertools.imap(int, c2_line)

        # ind_genotypes[ind_num] = itertools.imap(operator.add, c1_line, c2_line)
        ind_genotypes[ind_num] = itertools.izip(c1_line, c2_line)
        pass

    return itertools.izip(positions, itertools.izip(*ind_genotypes))


def process_ms_block_to_snp_list(ms_file, opts):

    """
    Process one ms block into a list of snp structures
    returns snps.
    """
    
    snps = ['hey'] * segsites

    for snp_num, (pos, genotypes) in enumerate(process_ms_block_to_genotypes(ms_file, opts.num_inds, opts.window_length)):

        snp_d = {}
        snp_d['pos'] = pos

        t_gt = list(itertools.imap(lambda j : genotypes[j], opts.updated_ind_indexed_to_orig_file_target))
        r_gt = list(itertools.imap(lambda j : genotypes[j], opts.updated_ind_indexed_to_orig_file_reference))
        genotypes = t_gt + r_gt

        snp_d['genotypes'] = genotypes
        snp_d['target'] = sum(t_gt) > 0
        snp_d['reference'] = sum(r_gt) > 0

        snps[snp_num] = snp_d
        
        pass


    return snps
    



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
