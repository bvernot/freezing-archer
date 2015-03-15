from __future__ import division
import sys, os, random, itertools
import fileinput
from operator import itemgetter
import argparse
#from numpy.random import binomial
#import tables
import re
import math


parser = argparse.ArgumentParser(description='Generate ms code for various models; optionally run the ms code with -r.')

parser.add_argument('-m', '--model', choices=['3p', '3pb', 'hg_null', 'hg_int_pyg', 'hg_int_san', 'hg_int_hdz', 'wall_int', 'wall_int_rev', 'asn_ea_aa', 'deme'], required = True)
parser.add_argument('-d', '--debug', action = 'store_true')
parser.add_argument('-dd', '--debug-2', action = 'store_true')
parser.add_argument('-ddd', '--debug-3', action = 'store_true')
parser.add_argument('-l', '--reglen', default=50000, type=int)
parser.add_argument('-p1', '--pop1-ninds', default=5, type=float)
parser.add_argument('-p2', '--pop2-ninds', default=5, type=float)
parser.add_argument('-p3', '--pop3-ninds', default=5, type=float)
parser.add_argument('-s1', '--pop1-size', default=None, type=float)
parser.add_argument('-s2', '--pop2-size', default=None, type=float)
parser.add_argument('-s2b', '--pop2-size-at-bottleneck', default=None, type=float)
parser.add_argument('-s3b', '--pop3-size-at-bottleneck', default=None, type=float)
parser.add_argument('-replace-ms', '--replace-ms-file', default=True)
parser.add_argument('-split', '--p1-p2-split', default=128000, type=int)
parser.add_argument('-split23', '--p2-p3-split', default=23000, type=float)
parser.add_argument('-struct', '--ancestral-structure', default=None, type=float, help='This is the percentage of each subpopulation (in the structured, ancestral population) that is made of migrants from the other population.')
parser.add_argument('-struct-join', '--struct-pop-join', default=None, type=int)
parser.add_argument('-is-struct', '--is-struct', default = 2, type = int, help='for ABC, sometimes we give struct params when we really don\'t want structure.  if that\'s the case, give the --is-struct 1 flag, and all struct params will be ignored.  It\'s really stupid, but I can\'t get ABCtoolbox to generate 1s and 0s, so instead 1 == False and 2 == True')
parser.add_argument('-N0', '--Nzero', default=10000, type=float)
parser.add_argument('-eur-bot', '--eur-bottleneck-size', default=91, type=float)
parser.add_argument('-gen-time', '--generation-time', default=20, type=int)
parser.add_argument('-arc-chrs', '--arc-chrs', default=1, type=int)
parser.add_argument('-eur-bot-time', '--eur-bottleneck-time', default=62000, type=float)
parser.add_argument('-target', '--target-pop', default=[2], type=int, nargs = '+')
parser.add_argument('-mig', '--migration-percent', default=[.01], type=float, nargs = '+')
parser.add_argument('-migt', '--migration-time', default=[1200*20], type=float, nargs = '+')
parser.add_argument('-migd', '--migration-duration', default=[20], type=float, nargs = '+')
parser.add_argument('-ast', '--arc-split-time', default=700000, type=float)
parser.add_argument('-recomb', '--recomb', default='1e-8')
parser.add_argument('-tbsfile', '--tbsfile', default=None)
parser.add_argument('-mu', '--mu', default=2.3e-8, type=float)
parser.add_argument('-mue8', '--mue8', default=None, type=float)
parser.add_argument('-n', '--num-samples', default=1, type=int)
parser.add_argument('-seg', '--seg-sites', default=None, type=int)
parser.add_argument('-wl', '--winlen', default=50000, type=int)
parser.add_argument('-ws', '--winstep', default=20000, type=int)
parser.add_argument('-iter', '--iteration', default=0, type=int)
parser.add_argument('-r', '--run', action='store_true')
parser.add_argument('-msage', '--use-msage-bin', action='store_true')
parser.add_argument('-ms', '--ms-bin', default = '/net/gs/vol2/home/bvernot/bin/ms')
parser.add_argument('-rc', '--run-and-calc', action='store_true')
parser.add_argument('-ra', '--arc-and-calc', action='store_true')
parser.add_argument('-T', '--report_trees', action='store_const', const='-T ', default='')
parser.add_argument('-o', '--output-file', type=argparse.FileType('w'), required = False, default = sys.stdout, help = 'output file (for ms output)')
#parser.add_argument('-pot', '--pipe-output-to', type=argparse.FileType('w'), required = False, default = sys.stdout, help = 'pipe all output to this file')
#parser.add_argument('-b', '--bottleneck', nargs = 4, , metavar=['pop', 'frac', 'start', 'end'], 'Add a population bottleneck to )
parser.add_argument('-seeds', '--seeds', default=None, type=int, nargs=3)

parser.add_argument('-tags', '--tags', required = False, default = [], nargs = '+')
parser.add_argument('-tagl', '--tag-labels', required = False, default = [], nargs = '+')

parser.add_argument('--deme-N', default=None, type=int)
parser.add_argument('--deme-arc-ea-dist', default=None, type=int)
parser.add_argument('--deme-size-ea', default=None, type=int)
parser.add_argument('--deme-ea-aa-dist', default=None, type=int)
parser.add_argument('--deme-size-aa', default=None, type=int)
parser.add_argument('--deme-mig', default=None, type=float)


args = parser.parse_args()

if len(args.tags) > 0:
    if len(args.tags) != len(args.tag_labels):
        print "error in tag logic!  must be the same length:"
        print args.tag_labels
        print args.tags
        sys.exit(-1)
        pass
    print ' '.join(['run_tag_labels'] + args.tag_labels)
    print ' '.join(['run_tags'] + args.tags)
    pass

# print args

num_tp = len(args.target_pop)
if num_tp != len(args.migration_percent) or num_tp != len(args.migration_time) or num_tp != len(args.migration_duration):
    print "target pop has to have the same length as all migration parameters:"
    print len(args.target_pop), len(args.migration_percent), len(args.migration_time), len(args.migration_duration)
    sys.exit(-1)
    pass

if args.use_msage_bin:
    args.ms_bin = '/net/gs/vol2/home/bvernot/bin/msage'
    pass

#sys.stdout = args.pipe_output_to

## set debug arguments
if args.debug_3: args.debug_2 = True
if args.debug_2: args.debug = True
if args.debug: print args

# for ABC, sometimes we give struct params when we really don't want structure.  if that's the case, give the --is-struct 0 flag, and all struct params will be ignored.
if args.is_struct == 1:
    args.ancestral_structure = None
    pass
if args.struct_pop_join == None:
    args.struct_pop_join = args.arc_split_time * 2
    pass

if args.mue8 != None:
    args.mu = args.mue8 / 100000000
    if args.debug: print 'setting mu!', args.mue8, args.mu
    pass

theta = 4 * args.Nzero * args.mu * args.reglen
if args.recomb != 'tbs':
    recomb = 4 * args.Nzero * float(args.recomb) * (args.reglen-1)
elif args.tbsfile == None:
    print 'cannot set recomb to tbs without -tbsfile!'
    sys.exit(-1)
else:
    recomb = 'tbs'
    pass


def scale_time(years, gen=args.generation_time):
    return years / gen / 4 / args.Nzero

def get_years(scaled_param, gen=args.generation_time):
    return scaled_param * gen * 4 * args.Nzero

def scale_migration(pct):
    return pct * 4 * args.Nzero

def scale_N(N):
    return N / args.Nzero

def get_alpha(Npresent, Npast, t):
    # t should be scaled time
    return -math.log(Npast / Npresent) / t

def get_N_at_t_units_ago(Npresent, alpha, t):
    # t should be scaled time
    return Npresent * math.exp(-alpha * t)


output_file_arg = ''
if args.output_file != sys.stdout:
    output_file_arg = '-o %s ' % (args.output_file.name)
    pass


if args.model == 'asn_ea_aa':

    # 1 is afr, 2 is eur, 3 is asn, 4 is introgressed (the commands here that are commented out are wenqing's commands for just afr and eur)

    # ./ms Indnum Repnum -I 2 Indnum1 Indnum2 -t 4Nu -r 4Nr Len 
    # -n 1 58.002735978 -n 2 70.041039672 
    # -eg 0 1 482.46 -eg 0 2 570.18 
    # -em 0 1 2 0.7310 -em 0 2 1 0.7310 
    # -eg 0.006997264 1 0 -eg 0.006997264 2 89.7668 
    # -en 0.006997264 1 1.98002736 -en 0.031463748 2 0.141176471 
    # -en 0.03146375 2 0.254582763 
    # -em 0.03146375 1 2 4.386 -em 0.03146375 2 1 4.386 
    # -em 0.069767442 1 2 0 -em 0.069767442 2 1 0 
    # -ej 0.069767442 2 1 -en 0.069767442 1 1.98002736 
    # -en 0.20246238 1 1
    

    if args.Nzero != 7310:
        print "ea_aa models require Nzero = 7310"
        sys.exit(-1)
        pass

    if args.generation_time != 25:
        print "ea_aa models require gen time = 25"
        sys.exit(-1)
        pass

    cmd = args.ms_bin + " %d %d " % (args.arc_chrs + args.pop1_ninds * 2 + args.pop2_ninds * 2 + args.pop3_ninds * 2, args.num_samples)

    if args.seg_sites != None:
        cmd += "-s %d " % args.seg_sites
    else:
        cmd += "-t %e " % theta
        pass
    cmd += "-r %s %d " % (recomb, args.reglen)
    if args.tbsfile != None:
        cmd = 'cat %s | %s' % (args.tbsfile, cmd)
        pass
    ## population sample sizes
    cmd += "-I 4 %d %d %d %d 0 " % (args.pop1_ninds * 2, args.pop2_ninds * 2, args.pop3_ninds * 2, args.arc_chrs)

    if args.debug:
        print get_years(0.006997264)
        print get_years(0.031463748)
        print get_years(0.03146375)
        print get_years(0.069767442)
        print get_years(0.20246238)
        print get_years(scale_time(args.arc_split_time))
        pass

    ## set archaic population size to 1500
    ## ARCHAIC POPULATION IS NOW 4
    cmd += "-n 4 %e " % scale_N(1500)

    ## for now, just use wenqing's parameters (should convert these into using my methods)
    cmd += '-n 1 58.002735978 -n 2 70.041039672 -n 3 187.55 '
    ## start of growth
    cmd += '-eg 0 1 482.46 -eg 0 2 570.18 -eg 0 3 720.23 '
    eur_N_at_start_of_second_growth = get_N_at_t_units_ago(70.041039672, 570.18, 0.006997264)
    asn_N_at_start_of_second_growth = get_N_at_t_units_ago(187.55, 720.23, 0.006997264)
    print "european scaled pop size at beginning of second growth:", eur_N_at_start_of_second_growth
    print "asian scaled pop size at beginning of second growth:", asn_N_at_start_of_second_growth
    ## set migration
    cmd += '-em 0 1 2 0.7310 -em 0 2 1 0.7310 ' # afr and eur
    cmd += '-em 0 1 3 0.228072 -em 0 3 1 0.228072 ' # afr and asn
    ######## MIGRATION BTWN EUR AND ASN: 3.11 * 0.2924 (3.11 is migration value from Gravel et al, and 0.2924 is the factor used to get the appropriate ms value)
    cmd += '-em 0 2 3 0.909364 -em 0 3 2 0.909364 ' # eur and asn

    ## continuing growth / ending growth (it's ok to have the -en after the -eg for pop 1, because we're setting the growth rate to 0 - actually, they're redundant)
    print "european scaled pop size at beginning of first growth:", get_N_at_t_units_ago(eur_N_at_start_of_second_growth, 89.7668,  0.031463748 - 0.006997264)
    print "asian scaled pop size at beginning of first growth:",    get_N_at_t_units_ago(asn_N_at_start_of_second_growth, 113.3896, 0.031463748 - 0.006997264)

    eur_alpha2 = get_alpha(eur_N_at_start_of_second_growth, args.pop2_size_at_bottleneck / args.Nzero, 0.031463748 - 0.006997264)
    print "european alpha necessary to go from Ne of", args.pop2_size_at_bottleneck, "to correct size at start of second growth:", eur_alpha2
    print "european scaled pop size at beginning of first growth (adjusted alpha):", get_N_at_t_units_ago(eur_N_at_start_of_second_growth, eur_alpha2, 0.031463748 - 0.006997264)

    asn_alpha2 = get_alpha(asn_N_at_start_of_second_growth, args.pop3_size_at_bottleneck / args.Nzero, 0.031463748 - 0.006997264)
    print "asian alpha necessary to go from Ne of", args.pop2_size_at_bottleneck, "to correct size at start of second growth:", asn_alpha2
    print "asian scaled pop size at beginning of first growth (adjusted alpha):", get_N_at_t_units_ago(asn_N_at_start_of_second_growth, asn_alpha2, 0.031463748 - 0.006997264)

    # cmd += '-eg 0.006997264 1 0 -eg 0.006997264 2 89.7668 -eg 0.006997264 3 113.3896 '
    cmd += '-eg 0.006997264 1 0 -eg 0.006997264 2 %e -eg 0.006997264 3 %e ' % (eur_alpha2, asn_alpha2)
    cmd += '-en 0.006997264 1 1.98002736 '

    ## bottleneck in europe and asia (replaced by join?)
    #cmd += '-en 0.031463748 2 0.141176471 '
    #cmd += '-en 0.031463748 3 0.07578659 ' # bottleneck is more severe in asia
    cmd += '-en 0.031463748 2 %e ' % scale_N(args.pop2_size_at_bottleneck)
    cmd += '-en 0.031463748 3 %e ' % scale_N(args.pop3_size_at_bottleneck) # bottleneck is more severe in asia?

    ## join eur and asns
    cmd += '-ej %e 3 2 ' % scale_time(args.p2_p3_split)
    # cmd += '-en %e 2 0.254582763 ' % scale_time(args.p2_p3_split)
    cmd += '-en %e 2 %e ' % (scale_time(args.p2_p3_split), scale_N(max(args.pop2_size_at_bottleneck, args.pop3_size_at_bottleneck)))

    ## reestablish migration between afr and eur+asn
    cmd += '-em %e 1 2 4.386 -em %e 2 1 4.386 ' % (scale_time(args.p2_p3_split), scale_time(args.p2_p3_split))

    ## end migration (unnecessary)
    # cmd += '-em 0.069767442 1 2 0 -em 0.069767442 2 1 0 '

    # and join into one population (51kya)
    cmd += '-ej %e 2 1 -en %e 1 1.98002736 ' % (scale_time(args.p1_p2_split), scale_time(args.p1_p2_split))

    # and collapse to ancestral population size
    cmd += '-en 0.20246238 1 1 '
    
    ## introgression for one generation, for specified percent, at specified time, to the specified population (from population 4)
    for mig_index in range(len(args.migration_percent)):
        if args.migration_percent[mig_index] > 0:
            sys.stderr.write('Simulating archaic introgression: into %d at %f rate at %f time for %f time.\n' % (args.target_pop[mig_index], 
                                                                                                                 args.migration_percent[mig_index], 
                                                                                                                 args.migration_time[mig_index], 
                                                                                                                 args.migration_duration[mig_index]))
            cmd += "-em %e %d 4 %e -em %e %d 4 0 " % (scale_time(args.migration_time[mig_index]), args.target_pop[mig_index], scale_migration(args.migration_percent[mig_index]), scale_time(args.migration_time[mig_index]+args.migration_duration[mig_index]), args.target_pop[mig_index])
        else:
            sys.stderr.write('Simulating no archaic introgression.\n')
            pass
        pass
    ## allow archaic split time to change
    cmd += "-ej %e 4 1 " % (scale_time(args.arc_split_time))
    
    ## output tree stuff
    cmd += "-L " + args.report_trees

    if args.seeds != None:
        cmd += "-seeds " + ' '.join([str(i) for i in args.seeds]) + ' '
        pass

    print "removed -ms-no-randomize from tp2 calculation (never added it to the tp3 calc, which was the only one where it was needed"
    arc_cmd = "python bin/get_pct_arc_per_ind_from_ms_file_new.py -f msfile -tp 2 3 -noref -regions -summary -prog 10 -iter %d" % args.iteration
#    arc_cmd = "python bin/get_pct_arc_per_ind_from_ms_file_new.py -f msfile -tp 2 3 -noref -hap-freq -summary -prog 10 -iter %d" % args.iteration

    # tmrca_cmd = "python /net/akey/vol1/home/bvernot/tishkoff/metrics/calculate_tmrca2.py region_s_star -ms - -ms-pops %d %d %d 1 -tp 2 " % (args.pop1_ninds * 2, 
    #                                                                                                                                         args.pop2_ninds * 2, 
    #                                                                                                                                         args.pop3_ninds * 2)
    # tmrca_cmd += "-ms-win 50000 50000 -ms-reglen 50000 -pops 1 2 -no-pearson -forward-ms-chr-sstar --forward-ms-file -no-sum "
    # tmrca_cmd += "| python /net/akey/vol1/home/bvernot/tishkoff/metrics/calculate_tmrca2.py region_s_star -ms - -ms-pops %d %d %d 1 -tp 3 " % (args.pop1_ninds * 2,
    #                                                                                                                                            args.pop2_ninds * 2, 
    #                                                                                                                                            args.pop3_ninds * 2)
    # tmrca_cmd += "-ms-win 50000 50000 -ms-reglen 50000 -pops 1 3 -no-pearson -forward-ms-chr-sstar --forward-ms-file -no-sum "
    # tmrca_cmd += " | " + arc_cmd


    if args.debug:
        print cmd
        # print tmrca_cmd
        pass

    pass


print >> sys.stderr, cmd.replace('/net/gs/vol2/home/bvernot/bin/ms', 'python ~/Dropbox/msVis/msVis.py') + ' -years -N0 %d -gen %d' % (args.Nzero, args.generation_time)



if args.run:
    from subprocess import call
    #print cmd
    #call(cmd.split())
    call(cmd, shell=True)

elif args.run_and_calc:
    from subprocess import call

    if args.output_file != sys.stdout:
        cmd += " > " + args.output_file.name + '.ms'
        print cmd

        if args.replace_ms_file: call(cmd, shell=True)
        print tmrca_cmd.replace('msfile', args.output_file.name + '.ms')
        call(tmrca_cmd.replace('msfile', args.output_file.name + '.ms'), shell=True)

    else:
        cmd += " | " + tmrca_cmd.replace('msfile', '-')
        print 'command:'
        print cmd
        call(cmd, shell=True)
        pass
    pass
elif args.arc_and_calc:
    from subprocess import call

    # if args.output_file != sys.stdout:
    #     cmd += " > " + args.output_file.name + '.ms'
    #     print cmd
    #     if args.replace_ms_file: call(cmd, shell=True)
    #     print arc_cmd.replace('msfile', args.output_file.name + '.ms')
    #     call(arc_cmd.replace('msfile', args.output_file.name + '.ms'), shell=True)
    # else:

    cmd += " | " + arc_cmd.replace('msfile', '-')
    # print cmd
    call(cmd, shell=True)

    # pass
    pass
else:
    print
    print "ms command:"
    print cmd

    print
    print "msVis command:"
    print cmd.replace('/net/gs/vol2/home/bvernot/bin/ms', 'python ~/Dropbox/msVis/msVis.py') + ' -years -N0 %d -gen %d' % (args.Nzero, args.generation_time)

    print
    print "use -r to run ms command"

    # print
    # print "tmrca command (set -tp appropriately):"
    # print tmrca_cmd
    pass
