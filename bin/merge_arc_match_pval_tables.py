from __future__ import division
from numpy import std
from collections import defaultdict, Counter
import arc_match_pval_tables
import cPickle
import itertools

local_debug = False


arc_match_table = cPickle.load(open('pval_match_tables.v1/n_match_table.chr%d.YRI.pickle' % 1))
# print "loading match tables.."
# arc_match_tables = []
# for c in range(1,23):
#     f = 'pval_match_tables.v1/n_match_table.chr%d.YRI.pickle' % c
#     print c, f
#     arc_match_tables.append(cPickle.load(open(f, 'rb')))
#     pass

# print "... finished loading match tables"


print_table = True
if print_table:
    print "bases_bin nsites sfs std nmatches count"
    for b in arc_match_table.keys():
        for p in arc_match_table[b].keys():
            for sfs in arc_match_table[b][p].keys():
                for sd in arc_match_table[b][p][sfs].keys():
                    for nm in arc_match_table[b][p][sfs][sd].keys():
                        print b,p,sfs,sd,nm,arc_match_table[b][p][sfs][sd][nm] #, nm/p
                        pass
                    pass
                pass
            pass
        pass
    pass


sys.exit()



print_lim = 2

for i in arc_match_table.keys()[:print_lim]:
    print i
    for j in arc_match_table[i].keys()[:print_lim]:
        print i,j
        for k in arc_match_table[i][j].keys()[:print_lim]:
            print i,j,k
            for l in arc_match_table[i][j][k].keys()[:print_lim]:
                print i,j,k,l,arc_match_table[i][j][k][l]
                pass
            pass
        pass
    pass


(my_bases_bin, my_sites, my_sfs, my_std_dev) = (32000, 31, 14, 97)

s_dist = [i for i in arc_match_table[my_bases_bin].keys() if my_sites-4 < i < my_sites+4]
print s_dist
null_hits = 0
total_hits = 0
for i in s_dist:
    hits = list(itertools.chain.from_iterable([[[arc_match_table[my_bases_bin][i][s1][s2]
                                                 for s2 in arc_match_table[my_bases_bin][i][s1].keys() if my_std_dev-10 < s2 < my_std_dev+10]
                                                for s1 in arc_match_table[my_bases_bin][i].keys() if my_sfs-10 < s1 < my_sfs+10]]))

    print hits
    pass



print 

new_table = arc_match_pval_tables.mydd3()

for b in set(itertools.chain(*(t.keys() for t in arc_match_tables))):
    print b
    for p in set(itertools.chain(*(t[b].keys() for t in arc_match_tables))):
        for sfs in set(itertools.chain(*(t[b][p].keys() for t in arc_match_tables))):
            for sd in set(itertools.chain(*(t[b][p][sfs].keys() for t in arc_match_tables))):
                # for i,t in enumerate(arc_match_tables):
                #     print 'table%d:' % i, b,p,sfs,sd,t[b][p][sfs][sd] #, nm/p
                #     pass
                new_table[b][p][sfs][sd] = reduce(lambda x,y : x+y, [t[b][p][sfs][sd] for t in arc_match_tables])
                # print '  merged', new_table[b][p][sfs][sd]
                pass
            pass
        pass
    pass
    
#                 for nm in set(itertools.chain(*(t[b][p][sfs][sd].keys() for t in arc_match_tables))):


cPickle.dump(new_table, open('tmp.pickle.chr1.merged', 'wb'), -1)

