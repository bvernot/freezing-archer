from __future__ import division
import sys, itertools, operator
import read_vcf

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
    
    snps = list()

    for snp_num, (pos, genotypes) in enumerate(process_ms_block_to_genotypes(ms_file, opts.num_inds, opts.window_length)):

        snp_d = {}
        snp_d['pos'] = pos

        t_gt = list(itertools.imap(lambda j : genotypes[j], opts.reference_individuals_indexed_to_orig_file))
        r_gt = list(itertools.imap(lambda j : genotypes[j], opts.target_individuals_indexed_to_orig_file))
        e_gt = list(itertools.imap(lambda j : genotypes[j], opts.exclude_individuals_indexed_to_orig_file))
        genotypes = t_gt + r_gt

        snp_d['genotypes'] = genotypes
        snp_d['target'] = sum(t_gt) > 0
        snp_d['reference'] = sum(r_gt) > 0

        snps.append(snp_d)
        
        pass


    return snps
    



#@profile
def ms_to_genotypes_windowed(ms_file, winlen, winstep, ms_ind_pop_file, opts, start = 0):

    ## process sample ID info / results are saved to opts
    read_ms_header(ms_file, opts)

    snps = []
    keep_reading = True
    chrom = 1
    winstart = start
    winend = winlen

    while keep_reading:

        snps = process_ms_block_to_snp_list(ms_file, opts)
        
        if snps == None: break
        
        yield (chrom, winstart, winend, snps)
        chrom += 1
        pass

    pass



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
    """
    Set necessary options - mostly the indices of target and reference populations.
    """

    opts.ms_command = ms_file.readline().strip().split()
    opts.ms_seeds = ms_file.readline().strip().split()

    num_chroms = int(opts.ms_command[1])
    num_pops = opts.ms_pop_sizes[0]

    ## I'M NOT SURE HOW ROBUST THIS IS!
    ####
    SHOULD WE REALLY SET ARC POPS INSTEAD OF CHRS?  THAT IS WHAT WE CARE MORE ABOUT.

    pop_list_chr = list(itertools.chain.from_iterable([i for _ in range(x)] for i,x in enumerate(opts.ms_pop_sizes[1:])))
    pop_mapping = {i:p for i,p in enumerate(pop_list_chr)}
    print pop_mapping
    ind_list = ['i%d' % i for i in range(opts.ms_num_diploid_inds)]

    ## check that number of chromosomes in ms file matches ms_pop_sizes
    if num_chroms != sum(opts.ms_pop_sizes[1:]):
        print "Mismatch between number of simulated chromosomes (%d) and population definitions given with --ms-pop-sizes ( %r == %d )" % (num_chroms, 
                                                                                                                                           opts.ms_pop_sizes[1:],
                                                                                                                                           sum(opts.ms_pop_sizes[1:]))
        sys.exit(-1)
        pass

    ## check that number of chromosomes in ms file matches 2*ms_num_diploid_inds + len(ms_archaic_chromosomes)
    if num_chroms != 2 * opts.ms_num_diploid_inds + len(opts.ms_archaic_chromosomes):
        print "Mismatch between number of simulated chromosomes (%d) and number of diploid / archaic individuals ( 2*%d + len(%r) )" % (num_chroms,
                                                                                                                                        opts.ms_num_diploid_inds,
                                                                                                                                        opts.ms_archaic_chromosomes)
        sys.exit(-1)
        pass

    ## check that a) the archaic chromosomes are at the end, and b) that all other populations have diploid individuals (even number of chrs)
    problematic_pops = [(pop_idx, pop_size) for pop_idx, pop_size in enumerate(opts.ms_pop_sizes[1:]) if pop_size % 2 != 0]
    print "problematic_pops", problematic_pops, opts.ms_archaic_chromosomes
    if problematic_pops:
        print "Mismatch between number of simulated chromosomes (%d) and number of diploid / archaic individuals ( 2*%d + len(%r) )" % (num_chroms,
                                                                                                                                        opts.ms_num_diploid_inds,
                                                                                                                                        opts.ms_archaic_chromosomes)
        sys.exit(-1)
        pass

    ## construct the individual to pop mapping, which usually is created from a file that looks like this:
    # sample  pop     super_pop       gender
    # HG00096 GBR     EUR     male
    # HG00097 GBR     EUR     female
    # HG00099 GBR     EUR     female
    # ind_pop_mapping = [l.strip().split()[:3] for l in f.readlines()]
    ## but for ms files we construct it out of the -I flag (here given by --ms-pop-sizes)

    ## I'M NOT SURE HOW ROBUST THIS IS - FOR ONE, IT REMOVES POPS WITH ONE CHR
    pop_list = list(itertools.chain.from_iterable([i for _ in range(x//2)] for i,x in enumerate(opts.ms_pop_sizes[1:])))

    print pop_list
    print ind_list
    print list(itertools.izip(pop_list, ind_list, ind_list))
    ind_pop_mapping = list(itertools.izip(pop_list, ind_list, ind_list))


    
    read_vcf.process_ind_pop_mapping(opts, ind_pop_mapping)

    
    ## could read trees here - do not currently do that

    line = ms_file.readline()
    while line != '' and not line.startswith('//'):
        line = ms_file.readline()
        pass


    opts.original_file_ind_ids = ['i%d' % i for i in range(opts.num_inds)]
    #opts.original_file_ind_index = {id:i for i,id in enumerate(opts.original_file_ind_ids)}
    return




# def read_ms_ind_pop_file(f, opts):

#     # sample  pop     super_pop       gender
#     # HG00096 GBR     EUR     male
#     # HG00097 GBR     EUR     female
#     # HG00099 GBR     EUR     female

#     header = f.readline()
#     opts.sample_to_pop = {}
#     opts.sample_to_superpop = {}
#     opts.all_pops = set()
#     opts.all_superpops = set()
#     opts.updated_reference_population_inds = set()
#     opts.updated_target_population_inds = set()

#     for line in f:
#         (sample, pop, superpop, _) = line.strip().split()
#         opts.sample_to_pop[sample] = pop
#         opts.sample_to_superpop[sample] = superpop
#         opts.all_pops.add(pop)
#         opts.all_superpops.add(superpop)

#         if pop      in opts.reference_populations: opts.updated_reference_population_inds.add(sample)
#         if superpop in opts.reference_populations: opts.updated_reference_population_inds.add(sample)
#         if pop      in opts.target_populations:    opts.updated_target_population_inds.add(sample)
#         if superpop in opts.target_populations:    opts.updated_target_population_inds.add(sample)

#         pass

#     opts.updated_target_population_inds.update(opts.target_individuals)
#     opts.updated_reference_population_inds.update(opts.reference_individuals)

#     opts.updated_ind_ids_target    = [i for i in opts.original_file_ind_ids if i in opts.updated_target_population_inds]
#     opts.updated_ind_ids_reference = [i for i in opts.original_file_ind_ids if i in opts.updated_reference_population_inds]
#     opts.updated_ind_indexed_to_orig_file_target    = [opts.original_file_ind_ids.index(ind) for ind in opts.updated_ind_ids_target]
#     opts.updated_ind_indexed_to_orig_file_reference = [opts.original_file_ind_ids.index(ind) for ind in opts.updated_ind_ids_reference]

#     opts.updated_ind_ids = opts.updated_ind_ids_target + opts.updated_ind_ids_reference
#     # opts.updated_ind_indexed_to_orig_file = opts.updated_ind_indexed_to_orig_file_target + opts.updated_ind_indexed_to_orig_file_reference

    
#     return# (sample_to_pop, sample_to_superpop)


