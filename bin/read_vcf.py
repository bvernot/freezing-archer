import sys, itertools


#@profile
def process_vcf_line_to_genotypes(line, opts):
    ## ignore comments
    ## header starts with just one #, comments start with ## - distinguish?
    ## we've already read the header (to get ind_ids)
    if line.lstrip().startswith('#'): return None
    
    snp_d = {}
    
    split_line = line.strip().split()
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096 HG00097
    snp_d['chrom'] = split_line[0]
    chrom = ('chr' if opts.vcf_has_illumina_chrnums else '') + snp_d['chrom']
    snp_d['pos'] = int(split_line[1])
    snp_d['ref'] = split_line[3]
    snp_d['alt'] = split_line[4]
    snp_d['qual'] = split_line[5]
    snp_d['filt'] = split_line[6]
    snp_d['info'] = split_line[7]
    format = split_line[8]

    ## remove non-biallelic sites
    if len(snp_d['ref']) != 1 or len(snp_d['alt']) != 1:
        if opts.debug: print "dropping snp because it is not biallelic, or is an indel", snp_d['chrom'], snp_d['pos'], snp_d['ref'], snp_d['alt']
        return None

    ## remove sites not in regions
    if opts.regions != None and not opts.regions.in_region_one_based(chrom, snp_d['pos']):
        if opts.debug: print "dropping snp because it is not in our regions definition:", snp_d['chrom'], snp_d['pos']
        return None
    if opts.regions != None and opts.debug: print "KEEPING SNP because it IS in our regions definition:", snp_d['chrom'], snp_d['pos']

    ## remove sites not in ancestral genome
    use_derived = opts.ancestral_bsg != None
    if use_derived and 'N' == opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos']):
        if opts.debug: print "dropping snp because it is not in the ancestral genome:", snp_d['chrom'], snp_d['pos'], opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos'])
        return None

    ## flip for derived?
    flip_for_derived = use_derived and snp_d['alt'].upper() == opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos'])
    if use_derived and opts.debug: print "checking derived from ancestral genome:", snp_d['chrom'], snp_d['pos'], \
            'ancestral =', opts.ancestral_bsg.get_base_one_based(chrom, snp_d['pos']), \
            'ref/alt =', snp_d['ref'], snp_d['alt'], 'flip_for_derived =', flip_for_derived, 'DERIVED' if flip_for_derived else ''
    

    ## get genotypes
    genotypes = split_line[9:]

    ## remove individual genotypes that aren't in the target or reference, and reorder (target first, then ref)
    ## (reordered ind ids in read_1kg_ind_pop_file)
    r_gt = list(itertools.imap(lambda j : genotypes[j], opts.reference_individuals_indexed_to_orig_file))
    t_gt = list(itertools.imap(lambda j : genotypes[j], opts.target_individuals_indexed_to_orig_file))
    e_gt = list(itertools.imap(lambda j : genotypes[j], opts.exclude_individuals_indexed_to_orig_file))


    if flip_for_derived:
        flip_map = {'.':'1|1', '0|0':'1|1', '1|0':'0|1', '0|1':'1|0', '1|1':'0|0'}
        t_gt_flipped = list(itertools.imap(lambda gt : flip_map[gt], t_gt))
        r_gt_flipped = list(itertools.imap(lambda gt : flip_map[gt], r_gt))
        e_gt_flipped = list(itertools.imap(lambda gt : flip_map[gt], e_gt))
        # print 'ORIGINAL GENOTYPES', t_gt, r_gt
        # print 'FLIPPED  GENOTYPES', t_gt_flipped, r_gt_flipped
        t_gt = t_gt_flipped
        r_gt = r_gt_flipped
        e_gt = e_gt_flipped

        ## SWAP REF AND ALT, so ref is "ancestral"
        # not sure if this is the best thing, because you might still want to know what the real ref/alt bases are.. but for now it works
        # (snp_d['ref'], snp_d['alt'])  = (snp_d['alt'], snp_d['ref'])
        pass

    genotypes = t_gt + r_gt
    t_gt = set(t_gt)
    r_gt = set(r_gt)

    e_gt = set(e_gt)
    if '0|1' in e_gt or '1|0' in e_gt or '1|1' in e_gt:
        if opts.debug: print "dropping snp because it's present in the EXCLUDE population"
        return None

    snp_d['target'] = '0|1' in t_gt or '1|0' in t_gt or '1|1' in t_gt
    # I DO NOT UNDERSTAND WHAT I WAS GETTING AT HERE...
    # snp_d['reference'] = not('0|0' in r_gt and len(r_gt) == 1)
    snp_d['reference'] = '0|1' in r_gt or '1|0' in r_gt or '1|1' in r_gt

    ## finally, remove variants that aren't in the target or reference
    if not snp_d['target'] and not snp_d['reference']:
        if opts.debug: print "dropping snp because it's NOT present in either the TARGET or REFERENCE"
        return None

    ## also remove variants that are fixed in target and reference
    if '1|1' in t_gt and len(t_gt) == 1 and '1|1' in r_gt and len(r_gt) == 1:
        if opts.debug: print "dropping snp because it's FIXED in both the TARGET and REFERENCE"
        return None


    ## get the location of genotypes, if necessary
    if True or ':' in format:
        gt_loc = format.split(':')
        if 'GT' not in gt_loc:
            print "bad format?  GT not found"
            print format
            print split_line
            sys.exit(-1)
            pass
        gt_loc = gt_loc.index('GT')
    
        ## strip out just genotypes
        # genotypes[:] = [gt.split(':')[gt_loc] for gt in genotypes]
        genotypes[:] = list(itertools.imap(lambda gt : gt.split(':')[gt_loc], genotypes))
        pass
    
    if opts.debug: print 'reading', snp_d['pos'], genotypes

    gt_map = {'0|0':0, '1|0':1, '0|1':1, '1|1':2}
    gt_map = {'.':0, '0|0':0, '1|0':1, '0|1':1, '1|1':2}
    hap_map1 = {'.':0, '0|0':0, '1|0':1, '0|1':0, '1|1':1}
    hap_map2 = {'.':0, '0|0':0, '1|0':0, '0|1':1, '1|1':1}
    gt_final = list(itertools.imap(lambda gt : gt_map[gt], genotypes))
    hap1_final = list(itertools.imap(lambda gt : hap_map1[gt], genotypes))
    hap2_final = list(itertools.imap(lambda gt : hap_map2[gt], genotypes))
    
    genotypes = gt_final
    if opts.debug: print genotypes
    snp_d['genotypes'] = genotypes
    snp_d['haplotypes_1'] = hap1_final
    snp_d['haplotypes_2'] = hap2_final

    snp_d['sfs_target'] = sum(genotypes[:opts.num_target])
    snp_d['sfs_reference'] = sum(genotypes[opts.num_target:opts.num_reference])
    if opts.neand_vcf != None:
        if flip_for_derived:
            # if the ref is derived (alt is ancestral), then check to see if neand has reference
            snp_d['arc_match'] = opts.neand_vcf.has_ref(snp_d['chrom'], snp_d['pos'])
        else:
            # otherwise, check to see if the alt matches neand
            snp_d['arc_match'] = snp_d['alt'].upper() in opts.neand_vcf.get_alts_one_based(snp_d['chrom'], snp_d['pos'])
            # print 'debug arc_match', snp_d['pos'], snp_d['alt'].upper(), opts.neand_vcf.get_alts_one_based(snp_d['chrom'], snp_d['pos'])
            pass
        snp_d['arc_genotype'] = opts.neand_vcf.get_genotype_one_based(snp_d['chrom'], snp_d['pos'], flip_for_derived)

        if opts.ancestral_bsg != None:
            snp_d['arc_is_derived'] = snp_d['arc_match']
            # if flip_for_derived:
            #     snp_d['arc_is_derived'] = opts.neand_vcf.has_ref(snp_d['chrom'], snp_d['pos'])
            # else:
            #     snp_d['arc_is_derived'] = snp_d['alt'].upper() in opts.neand_vcf.get_alts_one_based(snp_d['chrom'], snp_d['pos'])
            #     pass
            pass
            
        if opts.debug: print "neand match", snp_d['pos'], genotypes, snp_d['alt'].upper(), \
                opts.neand_vcf.get_ref_one_based(snp_d['chrom'], snp_d['pos']), \
                opts.neand_vcf.get_alts_one_based(snp_d['chrom'], snp_d['pos']), \
                opts.neand_vcf.get_bases_one_based(snp_d['chrom'], snp_d['pos']), \
                opts.neand_vcf.has_genotype_at_site(snp_d['chrom'], snp_d['pos']), \
                'DEBUG', \
                'ARCHAIC_MATCH' if snp_d['arc_match'] else '', \
                'DERIVED' if flip_for_derived else '', \
                'JACKPOT' if flip_for_derived and snp_d['arc_match'] and not opts.neand_vcf.has_ref(snp_d['chrom'], snp_d['pos']) else '', \
                'MINI_JACKPOT' if flip_for_derived and snp_d['arc_match'] else '', \
                'NEAND_IS_HET' if len(opts.neand_vcf.get_bases_one_based(snp_d['chrom'], snp_d['pos'])) == 2 else 'NEAND_IS_HOM'
                
        pass


    

    # ## do population munging?
    # present_pops = set([opts.sample_to_pop[opts.original_file_ind_ids[i]] for i in range(len(genotypes)) if genotypes[i] > 0])
    # present_superpops = set([opts.sample_to_superpop[opts.original_file_ind_ids[i]] for i in range(len(genotypes)) if genotypes[i] > 0])
    # present_inds = set([opts.original_file_ind_ids[i] for i in range(len(genotypes)) if genotypes[i] > 0])
    # if opts.debug: print present_pops, sum(genotypes)
    # if opts.debug: print present_superpops, sum(genotypes)
    
    # ## keep snps in the reference population..
    # ref_snp = False
    # for pop in opts.reference_populations:
    #     if pop in present_pops or pop in present_superpops:
    #         # if opts.debug: print "dropping snp because it IS present in the REFERENCE"
    #         ref_snp = True
    #         #return None
    #         break
    #     pass
        
    # ## keep snps in the reference individuals..
    # for ind in opts.reference_individuals:
    #     if ind in present_inds:
    #         # if opts.debug: print "dropping snp because it IS present in the REFERENCE"
    #         ref_snp = True
    #         #return None
    #         break
    #     pass
        
    # ## keep snps in the target population..
    # target_snp = False
    # for pop in opts.target_populations:
    #     if pop in present_pops or pop in present_superpops:
    #         target_snp = True
    #         break
    #     pass

    # ## also keep snps in the set of target individuals..
    # for ind in opts.reference_individuals:
    #     if ind in present_inds:
    #         target_snp = True
    #         break
    #     pass

    # ## finally, remove variants that aren't in the target or reference
    # if not (target_snp or ref_snp):
    #     if opts.debug: print "dropping snp because it's NOT present in either the TARGET or REFERENCE"
    #     return None

    # ## mark the populations that this snp is in (is this necessary?)
    # snp_d['populations'] = {}
    # for pop in opts.all_pops:
    #     snp_d['populations'][pop] = pop in present_pops
    #     pass
    # for pop in opts.all_superpops:
    #     snp_d['populations'][pop] = pop in present_superpops
    #     pass

    if opts.debug: print

    return snp_d


def vcf_to_genotypes(vcf_file):

    snps = []

    for line in vcf_file:
        
        snp_d = process_vcf_line_to_genotypes(line)
        if snp_d == None: continue
        snps.append(snp_d)

        pass

    return snps


#@profile
def vcf_to_genotypes_windowed(vcf_file, winlen, winstep, vcf_ind_pop_file, opts, start = 0):

    ## get sample ID order from the header (ignore comments) / results are saved to opts
    read_vcf_header(vcf_file, opts)

    ## process sample ID info / results are saved to opts
    read_1kg_ind_pop_file(vcf_ind_pop_file, opts)


    snps = []
    keep_reading = True
    chrom = None
    winstart = start
    winend = winlen

    while keep_reading:


        for line in vcf_file:

            if opts.debug: print 'line', line, 'line'
            
            snp_d = process_vcf_line_to_genotypes(line, opts)
            
            if snp_d == None: 
                if opts.debug: print
                continue

            chrom = snp_d['chrom']

            # if it's past the window, then yield the current set of snps
            # then adjust the window, prune snps, and check again if you should add it.. repeat..
            while snp_d['pos'] > winend:
                yield (chrom, winstart, winend, snps)
                winstart += winstep
                winend += winstep
                snps[:] = [s for s in snps if s['pos'] > winstart]
                pass
            
            # now it's finally definitely in the window!
            snps.append(snp_d)

        else:
            ## end of file
            keep_reading = False
            pass
        
        pass

    yield (chrom, winstart, winend, snps)



def process_ind_pop_mapping(opts, ind_pop_mapping):
    sample_to_pop = {}
    sample_to_superpop = {}

    all_pops = set()
    all_superpops = set()
    opts.reference_individuals = set(opts.reference_individuals)
    opts.target_individuals = set(opts.target_individuals)
    opts.exclude_individuals = set(opts.exclude_individuals)

    ## go through each line in the ind_pop file, and add the individual's ID (i.e., NA12078) to
    ##  opts.target_individuals, etc, if the pop or superpop matches
    pop_file_ind_id_order = []
    for sample, pop, superpop in ind_pop_mapping:

        sample_to_pop[sample] = pop
        sample_to_superpop[sample] = superpop
        all_pops.add(pop)
        all_superpops.add(superpop)
        pop_file_ind_id_order.append(sample)

        if pop      in opts.reference_populations: opts.reference_individuals.add(sample)
        if superpop in opts.reference_populations: opts.reference_individuals.add(sample)

        if pop      in opts.target_populations:    opts.target_individuals.add(sample)
        if superpop in opts.target_populations:    opts.target_individuals.add(sample)

        if pop      in opts.exclude_populations:   opts.exclude_individuals.add(sample)
        if superpop in opts.exclude_populations:   opts.exclude_individuals.add(sample)

        pass

    ## now sort them by order in the ind_pop file
    opts.target_individuals     = [i for i in pop_file_ind_id_order if i in opts.target_individuals]
    opts.exclude_individuals    = [i for i in pop_file_ind_id_order if i in opts.exclude_individuals]
    opts.reference_individuals  = [i for i in pop_file_ind_id_order if i in opts.reference_individuals]

    ## and make sure we have a mapping back to the order in the original vcf file
    opts.target_individuals_indexed_to_orig_file    = [opts.sample_index_in_original_file[ind] for ind in opts.target_individuals]
    opts.exclude_individuals_indexed_to_orig_file   = [opts.sample_index_in_original_file[ind] for ind in opts.exclude_individuals]
    opts.reference_individuals_indexed_to_orig_file = [opts.sample_index_in_original_file[ind] for ind in opts.reference_individuals]

    # print "opts.target_individuals", opts.target_individuals
    # print "opts.target_individuals", opts.target_individuals
    # print "pop_file_ind_id_order", pop_file_ind_id_order

    sample_ids = opts.target_individuals + opts.reference_individuals
    # opts.updated_ind_ids = opts.updated_ind_ids_target + opts.updated_ind_ids_reference
    # opts.updated_ind_indexed_to_orig_file = opts.updated_ind_indexed_to_orig_file_target + opts.updated_ind_indexed_to_orig_file_reference
    opts.num_target = len(opts.target_individuals)
    opts.num_reference = len(opts.reference_individuals)
    opts.num_samples = opts.num_target + opts.num_reference

    opts.target_indices = range(opts.num_target)
    opts.reference_indices = range(opts.num_target, opts.num_target+opts.num_reference)

    opts.get_id_from_sample_index = lambda ind: sample_ids[ind]
    opts.get_pop_from_sample_index = lambda ind: sample_to_pop[sample_ids[ind]]

    return


def read_1kg_ind_pop_file(f, opts):

    """
    Sets up important options, like the indices of reference and target individuals, ind to pop mappings, etc.
    Necessary to set:

    Indexes to the original VCF sample order:
    opts.reference_individuals_indexed_to_orig_file
    opts.target_individuals_indexed_to_orig_file
    opts.exclude_individuals_indexed_to_orig_file

    """

    # sample  pop     super_pop       gender
    # HG00096 GBR     EUR     male
    # HG00097 GBR     EUR     female
    # HG00099 GBR     EUR     female

    header = f.readline()

    ind_pop_mapping = [l.strip().split()[:3] for l in f.readlines()]

    process_ind_pop_mapping(opts, ind_pop_mapping)
    
    return# (sample_to_pop, sample_to_superpop)


def read_vcf_header(vcf_file, opts):

    line = vcf_file.readline()
    ## header starts with just one #, comments start with ##
    while line.lstrip().startswith('##'): line = vcf_file.readline()

    ##    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096 HG00097
    if not line.lstrip().startswith('#CHROM'):
        print "bad VCF header?"
        print line
        sys.exit(-1)
        pass

    # opts.original_file_ind_ids = line.strip().split()[9:]
    # opts.original_file_ind_index = {id:i for i,id in enumerate(opts.original_file_ind_ids)}
    opts.sample_index_in_original_file = {id:i for i,id in enumerate(line.strip().split()[9:])}
    return
