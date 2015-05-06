from __future__ import division
import sys, itertools, operator
import read_vcf
from custom_argparse import missingdict, vcf_class

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
def process_ms_block_to_genotypes(ms_file, num_inds, window_length, ms_chr, archaic_chromosomes_by_pop = []):

    """
    Process one ms block into positions and genotypes.
    Returns:
     - ((pos, gts), ..); a generator that returns consecutive positions and genotypes
     - a list of archaic vcf classes (from custom_argparse)
    """

    ##  assume that the // line has just been consumed


    ## get segsites
    line = ms_file.readline().strip()
    # print "READ LINE |%s|" % line
    
    while not line.startswith('segsites'):
        line = ms_file.readline()
        # print "READ LINE |%s|" % line
        if line == '':
            return None
        line = line.strip()
        pass

    segsites = int(line.split()[1])

    ## get position lists
    while not line.startswith('positions'):
        line = ms_file.readline().strip()
        pass
    
    # map the space [0,window_length)
    positions = list(itertools.imap(lambda s : int(float(s) * (window_length-1)), line.split()[1:]))

    ## read modern human genotypes, one ind at a time
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


    arc_vcfs = []
    for arc_pop_num, arc_pop_num_chrs in enumerate(archaic_chromosomes_by_pop):

        c1_line = ms_file.readline().strip()

        if arc_pop_num_chrs == 1:
            c2_line = c1_line
        elif arc_pop_num_chrs == 2:
            c2_line = ms_file.readline().strip()
        else:
            print "Currently can only handle one archaic individual per archaic population (i.e. 1 or 2 archaic chromosomes per population)"
            print "Archaic population %d has %d chromosomes." % (arc_pop_num, arc_pop_num_chrs)
            sys.exit(-1)
            pass
        
        c1_line = itertools.imap(int, c1_line)
        c2_line = itertools.imap(int, c2_line)

        vcf = {}
        # IF THIS IS CHANGED, ALSO CHANGE IT IN CUSTOM_ARGPARSE.PY! (should be set via a function..)
        vcf[ms_chr] = missingdict(lambda : (True, '0/0', 'N', ['N']))

        for pos, (g1, g2) in itertools.izip(positions, itertools.izip(c1_line, c2_line)):
            
            # IF THIS IS CHANGED, ALSO CHANGE IT IN CUSTOM_ARGPARSE.PY! (should be set via a function..)
            # (has_ref, genotype ('0/0', etc), ref_base, alt_bases)
            vcf[ms_chr][pos] = (g1 == 0 or g2 == 0, '%d/%d' % (g1, g2), '0', '1')
            pass
        
        arc_vcfs.append(vcf)
        pass

    return [itertools.ifilter(lambda x: sum([item for sublist in x[1] for item in sublist]) > 0, itertools.izip(positions, itertools.izip(*ind_genotypes))), arc_vcfs]


def process_ms_block_to_snp_list(ms_file, opts):
    
    """
    Process one ms block into a list of snp structures
    returns snps.
    """
    
    if opts.debug: print "PROCESSING MS BLOCK %d" % opts.current_ms_chrom_index 
    
    snps = list()
    opts.current_ms_chrom_index += 1
    opts.current_ms_chrom = 'ms%d' % opts.current_ms_chrom_index

    arc_pop_sizes = [len(c) for c in opts.ms_archaic_chromosomes_by_pop]
    genotype_results = process_ms_block_to_genotypes(ms_file,
                                                     opts.ms_num_diploid_inds,
                                                     opts.ms_simulated_region_length,
                                                     opts.current_ms_chrom,
                                                     arc_pop_sizes)

    if genotype_results == None:
        return None

    ms_snp_iterator, archaic_vcfs = genotype_results
    
    if len(archaic_vcfs) > 0:
        if len(archaic_vcfs) > 1: sys.stderr.write( "WARNING, ONLY TAKING ONE ARCHAIC AT THIS POINT\n" )
        opts.archaic_vcf = vcf_class(archaic_vcfs[0])
    else:
        opts.archaic_vcf = vcf_class({})
        opts.archaic_vcf.init_chrom(opts.current_ms_chrom)
        pass

    # print opts.target_individuals_indexed_to_orig_file
    
    for snp_num, (pos, genotypes) in enumerate(ms_snp_iterator):

        snp_d = {}
        snp_d['chrom'] = opts.current_ms_chrom
        snp_d['pos'] = pos

        t_gt = list(itertools.imap(lambda j : genotypes[j], opts.target_individuals_indexed_to_orig_file))
        r_gt = list(itertools.imap(lambda j : genotypes[j], opts.reference_individuals_indexed_to_orig_file))
        e_gt = list(itertools.imap(lambda j : genotypes[j], opts.exclude_individuals_indexed_to_orig_file))
        genotypes = t_gt + r_gt

        # print t_gt
        # print r_gt

        snp_d['ref'] = '0'
        snp_d['alt'] = '1'
        snp_d['genotypes'] = genotypes
        snp_d['target'] = (1,1) in t_gt or (0,1) in t_gt or (1,0) in t_gt
        snp_d['reference'] = (1,1) in r_gt or (0,1) in r_gt or (1,0) in r_gt
        snp_d['haplotypes_1'] = [i for i,_ in genotypes]
        snp_d['haplotypes_2'] = [j for _,j in genotypes]
        snp_d['arc_match'] = snp_d['alt'] in opts.archaic_vcf.get_alts_one_based(snp_d['chrom'], snp_d['pos'])

        # turn into derived count, not haplotypes
        snp_d['genotypes'] = [sum(gt) for gt in genotypes]

        snp_d['sfs_target'] = sum(snp_d['genotypes'][:opts.num_target])
        snp_d['sfs_reference'] = sum(snp_d['genotypes'][opts.num_target:opts.num_reference])
        
        snps.append(snp_d)
        
        pass

    # print snps
    return snps




#@profile
def ms_to_genotypes_windowed(ms_file, winlen, winstep, ms_ind_pop_file, opts, start = 0):

    ## process sample ID info / results are saved to opts
    read_ms_header(ms_file, opts)

    snps = []
    keep_reading = True
    winstart = start
    winend = winlen
    opts.current_ms_chrom_index = 0

    while keep_reading:

        snps = process_ms_block_to_snp_list(ms_file, opts)
        
        if snps == None: break

        for strt in range(0, opts.ms_simulated_region_length - winlen + 1, winstep):
            end = strt + winlen
            yield (opts.current_ms_chrom, strt, end, [s for s in snps if s['pos'] >= strt and s['pos'] < end])
            pass

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

    ## check that the ms-pop-sizes argument is well-formatted
    if opts.vcf_is_ms_file and opts.ms_pop_sizes != None and opts.ms_pop_sizes[0] != len(opts.ms_pop_sizes) - 1:
        print "--ms-pop-sizes is not well-formatted (%d pops required, but %d given).  Should follow -I from ms." % (opts.ms_pop_sizes[0], 
                                                                                                                     len(opts.ms_pop_sizes) - 1)
        sys.exit(-1)
        pass

    ## I'M NOT SURE HOW ROBUST THIS IS!
    ####
    # SHOULD WE REALLY SET ARC POPS INSTEAD OF CHRS?  THAT IS WHAT WE CARE MORE ABOUT.

    pop_list_chr = list(itertools.chain.from_iterable([i for _ in range(x)] for i,x in enumerate(opts.ms_pop_sizes[1:])))
    pop_mapping = {i:p for i,p in enumerate(pop_list_chr)}
    # print pop_list_chr
    # print pop_mapping
    ind_list = ['i%d' % i for i in range(opts.ms_num_diploid_inds)]

    ## check that the archaic population(s) are in the population list
    for arc_pop in opts.ms_archaic_populations:
        if arc_pop not in pop_list_chr:
            print "Archaic population %d is not valid!  Must be in: %r (determined from --ms-pop-sizes)" % (arc_pop, set(pop_list_chr))
            sys.exit(-1)
            pass
        pass

    ## check that number of chromosomes in ms file matches ms_pop_sizes
    if num_chroms != sum(opts.ms_pop_sizes[1:]):
        print ("Mismatch between number of simulated chromosomes (%d) and" + \
                   "population definitions given with --ms-pop-sizes ( %r == %d )") % (num_chroms, 
                                                                                       opts.ms_pop_sizes[1:],
                                                                                       sum(opts.ms_pop_sizes[1:]))
        sys.exit(-1)
        pass
    
    
    ####
    ## set up archaic chromosomes, both flat and by population

    opts.ms_archaic_chromosomes_by_pop = [[i for i,p in enumerate(pop_list_chr) if p == pop] for pop in opts.ms_archaic_populations]
    # print opts.ms_archaic_chromosomes_by_pop
    opts.ms_archaic_chromosomes_flat = [item for sublist in opts.ms_archaic_chromosomes_by_pop for item in sublist]
    # print opts.ms_archaic_chromosomes_flat

    ## check that number of chromosomes in ms file matches 2*ms_num_diploid_inds + len(opts.ms_archaic_chromosomes)
    if num_chroms != 2 * opts.ms_num_diploid_inds + len(opts.ms_archaic_chromosomes_flat):
        print "Mismatch between number of simulated chromosomes (%d) and number of diploid / archaic individuals ( 2*%d + len(%r) )" % (num_chroms,
                                                                                                                                        opts.ms_num_diploid_inds,
                                                                                                                                        opts.ms_archaic_chromosomes_flat)
        sys.exit(-1)
        pass

    ## check that a) the archaic chromosomes are at the end, and b) that all other populations have diploid individuals (even number of chrs)
    pop_indices = range(opts.ms_pop_sizes[0])
    if len(opts.ms_archaic_populations) > 0 and sorted(pop_indices[-len(opts.ms_archaic_populations):]) != sorted(opts.ms_archaic_populations):
        print "Archaic populations must be the last simulated populations (may no longer actually be important..)."
        print "Archaic pops:", opts.ms_archaic_populations
        print "Specified pops:", pop_indices
        print sorted(pop_indices[-len(opts.ms_archaic_populations):])
        print sorted(opts.ms_archaic_populations)
        sys.exit(-1)
        pass

    # print pop_indices[:len(opts.ms_archaic_populations)]
    
    problematic_pops = [(pop_idx, pop_size) for pop_idx, pop_size in enumerate(opts.ms_pop_sizes[1:-len(opts.ms_archaic_populations)]) if pop_size % 2 != 0]
    # print "problematic_pops", problematic_pops

    if len(problematic_pops) != 0:
        print "Non-archaic populations must be diploid.  Number of chromosomes per population: %r" % opts.ms_pop_sizes[1:-len(opts.ms_archaic_populations)]
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
    pop_list = list(itertools.chain.from_iterable([str(i) for _ in range(x//2)] for i,x in enumerate(opts.ms_pop_sizes[1:])))

    # print pop_list
    # print ind_list
    # print list(itertools.izip(ind_list, pop_list, pop_list))
    ind_pop_mapping = list(itertools.izip(ind_list, pop_list, pop_list))


    # save sample index (easy for ms files, but needed for read_vcf.process_ind_pop_mapping)
    opts.sample_index_in_original_file = {('i%d' % i):i for i in range(opts.ms_num_diploid_inds)}


    ## get correct individual order, etc
    read_vcf.process_ind_pop_mapping(opts, ind_pop_mapping)
    
    
    ## could read trees here - do not currently do that

    ## eat lines until you get to a "//", which indicates the start of an ms block
    line = ms_file.readline()
    # print "READ LINE |%s|" % line
    while line != '' and not line.startswith('//'):
        line = ms_file.readline()
        # print "READ LINE |%s|" % line
        pass


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


