from __future__ import division
import sys, itertools, operator
import read_vcf
from custom_argparse import missingdict
import read_archaic_vcf
from test_parse_tree import read_tree4, find_node, tree_to_str, get_terminals, get_path_to_root, get_dist_btwn_nodes
from collections import defaultdict



def extract_introgressed_chrs_from_tree(ms_tree, arc_chrs, join_time, print_join_time = False):

    debug = False
    
    if debug: print ms_tree
    # get bases including bracket, i.e., [142]
    bases = ms_tree[:ms_tree.index(']')+1]
    # print bases, ms_tree
    
    # trim the bases from the tree (using the length of bases+brackets)
    ms_tree = ms_tree[len(bases):]
    
    # get the int from within the brackets
    bases = int(bases[1:-1])

    # arbitrarily look for the first archaic chromosome
    arc_chr = arc_chrs[0]
    split_tree = ms_tree.split('%d:' % arc_chr)
    if len(split_tree) != 2:
        print 'error in finding branch len for archaic chr %d' % arc_chr
        print split_tree
        print ms_tree
        sys.exit(-1)
        pass
    #bls2 = re.findall(r'[\d\.]+', split_tree[1])
    end_dist = min((x if x >=0 else 10000000000 for x in (split_tree[1].find(','), split_tree[1].find('('), split_tree[1].find(')'))))
    bls = split_tree[1][:end_dist]
    #print 'bls', bls, bls2[0], end_dist, split_tree[1][:10], (split_tree[1].find(','), split_tree[1].find('('), split_tree[1].find(')')), min((split_tree[1].find(','), split_tree[1].find('('), split_tree[1].find(')')))
    bl = float(bls)
    ## if the branch leading to the archaic chromosome is longer than the join time of the 
    ##   archaic and modern populations, then there's no way to identify introgressed sequence
    if debug: print 'archaic node branch length vs speciation time', bl, join_time

    if print_join_time: print "JOIN", bl

    if bl >= join_time:
        ## originally was a continue
        return (bases, [])
    
    tree = read_tree4(ms_tree)
    #print tree
    #print tree_to_str(tree)

    arc_node = find_node(tree, str(arc_chr))
    #print arc_node, arc_chr
    #print [n['name'] for n in get_path_to_root(arc_node)]
    #print [get_dist_btwn_nodes(n, arc_node) for n in get_path_to_root(arc_node)]
    # this returns nodes in decending order - so as soon as we get a node (c) below the join time, we can just return there
    for c in get_path_to_root(arc_node):
        dist = get_dist_btwn_nodes(c, arc_node)
        # if arc_node != c: print c, dist, join_time
        if debug: print dist, [int(n['name']) for n in get_terminals(c)]
        if arc_node != c and dist < join_time:
            terminal_chrs = get_terminals(c)
            terminal_chrs = [int(n['name']) for n in terminal_chrs]
            if debug: print 'introgressed (or archaic) chromosomes:', len(terminal_chrs), terminal_chrs
            return (bases, [n for n in terminal_chrs if n not in arc_chrs])
        pass
    return (bases, [])


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
def process_ms_block_to_genotypes(ms_file, num_inds, window_length, ms_chr, archaic_chromosomes_by_pop = [], ms_archaic_populations_join_times = []):

    """
    Process one ms block into positions and genotypes.
    Returns:
     - ((pos, gts), ..); a generator that returns consecutive positions and genotypes
     - a list of archaic vcf classes (from custom_argparse)
    """

    ##  assume that the // line has just been consumed

    ## read trees here

    line = ms_file.readline().strip()
    # print "READ LINE PROCESS GENOTYPES - trees? |%s|" % line

    chr_intr_regions = defaultdict(list)
    current_pos = 0
    while line.startswith('['):
        # this currently expects a list of chromosomes to look for, but *it only looks for the first one*
        # so we send just the first population worth of chromosomes
        # also, we are dealing with 0 based chrs, which is not how the trees in ms are formatted, so be careful!
        # I'm adding 1 to the archaic chromosomes, to look for them in the trees
        # and subtracting one from the intr_chrs, to make them 0 based
        (bases, intr_chrs) = extract_introgressed_chrs_from_tree(line, [c+1 for c in archaic_chromosomes_by_pop[0]], ms_archaic_populations_join_times[0])
        intr_chrs = [c-1 for c in intr_chrs]
        # print bases, intr_chrs
        
        # we might not have any introgressed chromosomes in this particular region
        num_intr_chrs = len(intr_chrs)
        if num_intr_chrs > 0:
            ## the regions introgressed for each chr (chr_intr_regions[chrom])
            for chrom in intr_chrs:
                if len(chr_intr_regions[chrom]) > 0 and chr_intr_regions[chrom][-1][1] == current_pos:
                    chr_intr_regions[chrom][-1][1] = current_pos+bases
                else:
                    chr_intr_regions[chrom].append([current_pos, current_pos+bases])
                    pass
                # print 'adding introgression information about chromosome', chrom, 'for region', (current_pos, current_pos+bases), chr_intr_regions[chrom]
                pass
            pass
        current_pos += bases

        line = ms_file.readline().strip()
        # print "READ LINE PROCESS GENOTYPES - trees? |%s|" % line
        pass

    # print "DONE WITH TREES"


    # for chrom in chr_intr_regions:
    #     print 'INTR', chrom, chr_intr_regions[chrom]
    #     pass

    ## get segsites
    
    while not line.startswith('segsites'):
        line = ms_file.readline()
        # print "READ LINE PROCESS GENOTYPES |%s|" % line
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
    arc_pop_sizes = [len(c) for c in archaic_chromosomes_by_pop]
    for arc_pop_num, arc_pop_num_chrs in enumerate(arc_pop_sizes):

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

        vcf = read_archaic_vcf.vcf_class()
        # IF THIS IS CHANGED, ALSO CHANGE IT IN CUSTOM_ARGPARSE.PY! (should be set via a function..)
        vcf.init_chrom(ms_chr)

        for pos, (g1, g2) in itertools.izip(positions, itertools.izip(c1_line, c2_line)):
            
            # IF THIS IS CHANGED, ALSO CHANGE IT IN CUSTOM_ARGPARSE.PY! (should be set via a function..)
            # (chr, pos, genotype ('0/0', etc), ref_base, alt_base, anc_base)
            # vcf.add_site(ms_chr, pos, gt, ref, alt, anc)
            vcf.add_site(ms_chr, pos, '%d/%d' % (g1, g2), '0', '1', '0')
            pass
        
        arc_vcfs.append(vcf)
        pass

    ## now consume until we get to EOF or //
    while line != '' and not line.startswith('//'):
        line = ms_file.readline()
        # print "READ LINE |%s|" % line
        pass


    return [itertools.ifilter(lambda x: sum([item for sublist in x[1] for item in sublist]) > 0, 
                              itertools.izip(positions, itertools.izip(*ind_genotypes))), 
            arc_vcfs,
            chr_intr_regions]


def process_ms_block_to_snp_list(ms_file, opts):
    
    """
    Process one ms block into a list of snp structures
    returns snps.
    """
    
    if opts.debug: print "PROCESSING MS BLOCK %d" % opts.current_ms_chrom_index 
    
    snps = list()
    opts.current_ms_chrom_index += 1
    opts.current_ms_chrom = 'ms%d' % opts.current_ms_chrom_index

    genotype_results = process_ms_block_to_genotypes(ms_file,
                                                     opts.ms_num_diploid_inds,
                                                     opts.ms_simulated_region_length,
                                                     opts.current_ms_chrom,
                                                     opts.ms_archaic_chromosomes_by_pop,
                                                     opts.ms_archaic_populations_join_times)

    if genotype_results == None:
        return None

    ms_snp_iterator, archaic_vcfs, chr_intr_regions = genotype_results

    # print 'archaicVCF', archaic_vcfs
    
    if len(archaic_vcfs) > 0:
        if len(archaic_vcfs) > 1: sys.stderr.write( "WARNING, ONLY TAKING ONE ARCHAIC AT THIS POINT\n" )
        opts.archaic_vcf = archaic_vcfs[0]
    else:
        opts.archaic_vcf = read_archaic_vcf.vcf_class()
        opts.archaic_vcf.init_chrom(opts.current_ms_chrom)
        pass

    # print opts.target_individuals_indexed_to_orig_file
    
    # print 'tar', opts.target_individuals_indexed_to_orig_file
    # print 'ref', opts.reference_individuals_indexed_to_orig_file
    # print 'excl', opts.exclude_individuals_indexed_to_orig_file

    ## before remapping introgressed regions, potentially output all introgressed regions in bed format
    for intr_chrom in sorted(chr_intr_regions):
        for s,e in chr_intr_regions[intr_chrom]:
            # get the ind: intr_chrom // 2
            # and the hap: intr_chrom % 2
            print '\t'.join(str(s) for s in [opts.current_ms_chrom, s, e, 
                                             '%s_i%d' % (opts.current_ms_chrom, intr_chrom // 2), 
                                             '%s_i%d_%d' % (opts.current_ms_chrom, intr_chrom // 2, intr_chrom % 2), 
                                             intr_chrom, 
                                             'INTROGRESSED_HAP'])
            pass
        pass
    
    ## have to remap chr_intr_regions to inds and haps (using opts.target_individuals_indexed_to_orig_file, etc)
    chr_intr_regions_remap = {}
    # print opts.target_individuals_indexed_to_orig_file + opts.reference_individuals_indexed_to_orig_file
    for i,o in enumerate(opts.target_individuals_indexed_to_orig_file + opts.reference_individuals_indexed_to_orig_file):
        # print "setting ind %d with original chr positions %d,%d" % (i, o*2, o*2+1)
        # print chr_intr_regions[o*2], chr_intr_regions[o*2+1]
        chr_intr_regions_remap[i] = {}
        chr_intr_regions_remap[i][1] = chr_intr_regions[o*2]
        chr_intr_regions_remap[i][2] = chr_intr_regions[o*2+1]
        pass

    # print 'orig intr regions:', chr_intr_regions
    # for ind in chr_intr_regions_remap:
    #     print 'remap intr:', ind, chr_intr_regions_remap[ind]
    #     pass
    
    is_introgressed_fn = lambda ind, hap, pos : sum(s <= pos <= e for s,e in chr_intr_regions_remap[ind][hap]) > 0
    # print "is_introgressed(8, 1000)", is_introgressed(8, 1000)
    # print "is_introgressed(4, 3856)", is_introgressed(4, 3856)
    # print "is_introgressed(4, 3855)", is_introgressed(4, 3855)
    
    
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
        snp_d['haplotypes_1_intr'] = [bool(i) and is_introgressed_fn(ind, 1, pos) for ind,(i,_) in enumerate(genotypes)]
        snp_d['haplotypes_2_intr'] = [bool(j) and is_introgressed_fn(ind, 2, pos) for ind,(_,j) in enumerate(genotypes)]
        snp_d['arc_match'] = opts.archaic_vcf.has_derived(snp_d['chrom'], snp_d['pos'])

        # turn into derived count, not haplotypes
        snp_d['genotypes'] = [sum(gt) for gt in genotypes]

        snp_d['sfs_target'] = sum(snp_d['genotypes'][i] for i in opts.target_indices)
        snp_d['sfs_reference'] = sum(snp_d['genotypes'][i] for i in opts.reference_indices)
        
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

    ## these are 0 based chrs, which is not how the trees in ms are formatted, so be careful!
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

    ## these are 0 based chrs, which is not how the trees in ms are formatted, so be careful!
    opts.ms_archaic_chromosomes_by_pop = [[i for i,p in enumerate(pop_list_chr) if p == pop] for pop in opts.ms_archaic_populations]
    # print " opts.ms_archaic_chromosomes_by_pop",  opts.ms_archaic_chromosomes_by_pop
    ## these are 0 based chrs, which is not how the trees in ms are formatted, so be careful!
    opts.ms_archaic_chromosomes_flat = [item for sublist in opts.ms_archaic_chromosomes_by_pop for item in sublist]
    # print " opts.ms_archaic_chromosomes_flat", opts.ms_archaic_chromosomes_flat

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

    # ####
    # ## parse ms command line to get the archaic join time [CURRENTLY ONLY WORKS WITH ONE ARCHAIC POP]
    # if len(opts.ms_archaic_populations) > 0:
    #     arc_pop = opts.ms_archaic_populations[0]
    #     print "looking for arc pop", arc_pop
    #     join_time = 0
    #     join_time_count = 0
    #     mig_pop_count = 0
    #     introgression_times = defaultdict(list)
        
    #     for i, flag in enumerate(opts.ms_command):
    #         if opts.debug and flag == '-ej': print flag, str(opts.arc_pop), ms_params[i+2:i+4]
    #         if flag == '-ej' and (str(arc_pop) == ms_params[i+2] or str(arc_pop) == ms_params[i+3]):
    #             join_time = float(ms_params[i+1])
    #             join_time_count += 1
    #             # print "Archaic population join time:", join_time
    #             pass
    #         if flag == '-em' and str(opts.arc_pop_ms) == ms_params[i+3]:
    #             parsed_target_pop = int(ms_params[i+2])
    #             introgression_times[parsed_target_pop].append(float(ms_params[i+1]))
    #             # print introgression_times
    #             mig_pop_count += 1
    #             pass
    #         pass

    #     if opts.debug: print join_time, join_time_count, parsed_target_pop, mig_pop_count
    #     if join_time_count != 1:
    #         print 'Error in parsing join time!'
    #         print join_time, join_time_count
    #         sys.exit(-1)
    #         pass
    #     pass
        
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


