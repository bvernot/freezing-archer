from __future__ import division
import sys, os, random, itertools, tempfile, cStringIO
import fileinput, subprocess
from operator import itemgetter
import argparse
#from numpy.random import binomial
#import tables
import re
import math
#from Bio import Phylo
from StringIO import StringIO
from collections import Counter, defaultdict





def read_tree(tree, debug = False):
    tree = tree.strip()
    if tree.endswith(';'): tree = tree[:-1]
    node = {'name':'', 'leaf':False, 'parent':None, 'dist':0}

    ## if this is a leaf node
    if not '(' in tree:
        tree = tree.split(':')
        node['name'] = tree[0]
        node['leaf'] = True
        if len(tree) > 1: node['dist'] = float(tree[1])
        return node

    ## if this node has a name
    if tree[0] != '(':
        node['name'] = tree[0:tree.index('(')]
        pass

    ## if this node has a dist
    if tree[-1] != ')':
        node['dist'] = float(tree[tree.rindex('):')+2:])
        pass

    tree = tree[tree.index('(')+1:tree.rindex(')')]
    splittree = tree.split(',')
    if debug: print splittree

    paren_counts = [s.count('(') - s.count(')') for s in splittree]

    for comma in range(1,len(paren_counts)):
        if debug: print 'considering comma', comma
        if sum(paren_counts[:comma]) == 0 and sum(paren_counts[comma:]) == 0:
            node['left'] = read_tree(','.join(splittree[:comma]))
            node['right'] = read_tree(','.join(splittree[comma:]))
            pass
        pass

    node['left']['parent'] = node
    node['right']['parent'] = node
    
    return node


# def read_tree2_old(tree, debug = False, splittree = None, paren_counts = None):
    
#     if splittree == None:
        
#         tree = tree.strip()
#         if tree.endswith(';'): tree = tree[:-1]
        
#         splittree = [s.strip() for s in tree.split(',')]
#         paren_counts = [s.count('(') - s.count(')') for s in splittree]
        
#         if debug: print splittree
#         pass
    
#     node = {'name' : '',
#             'leaf' : False,
#             'parent' : None,
#             'dist' : 0}
    
#     first = splittree[0]
#     last = splittree[-1]
    
#     ## if this is a leaf node
#     if len(splittree) == 1:
#         tree = first.split(':')
#         node['name'] = tree[0]
#         node['leaf'] = True
#         if len(tree) > 1: node['dist'] = float(tree[1])
#         return node
    
#     ## if this node has a name
#     if first[0] != '(':
#         node['name'] = first[0:first.index('(')]
#         pass
    
#     ## if this node has a dist
#     if last[-1] != ')':
#         node['dist'] = float(last[last.rindex('):')+2:])
#         pass
    
#     ## remove outer parens
#     splittree[0] = first[first.index('(')+1:]
#     splittree[-1] = last[:last.rindex(')')]
#     paren_counts[0] -= 1
#     paren_counts[-1] += 1

#     for comma in range(1,len(paren_counts)):
#         if debug: print 'considering comma', comma
#         if sum(paren_counts[:comma]) == 0 and sum(paren_counts[comma:]) == 0:
#             node['left'] = read_tree2('', splittree = splittree[:comma], paren_counts = paren_counts[:comma])
#             node['right'] = read_tree2('', splittree = splittree[comma:], paren_counts = paren_counts[comma:])
#             break
#             pass
#         pass
    
#     node['left']['parent'] = node
#     node['right']['parent'] = node
    
#     return node


def read_tree2(tree, debug = False, splittree = None, paren_counts = None, paren_sums = None):
    
    if splittree == None:
        
        tree = tree.strip()
        if tree.endswith(';'): tree = tree[:-1]
        
        splittree = [s.strip() for s in tree.split(',')]
        paren_counts = [s.count('(') - s.count(')') for s in splittree]
        # paren_sums = [[sum(paren_counts[:comma]), sum(paren_counts[comma:])] for comma in range(1,len(paren_counts))]
        
        if debug: print splittree
        pass
    
    node = {'name' : '',
            'leaf' : False,
            'parent' : None,
            'dist' : 0}
    
    first = splittree[0]
    last = splittree[-1]

    # print splittree, first, last

    def decorate(node, first, last, splittree):
    ## if this is a leaf node
        if len(splittree) == 1:
            tree = first.split(':')
            node['name'] = tree[0]
            node['leaf'] = True
            if len(tree) > 1: node['dist'] = float(tree[1])
            return
        
    ## if this node has a name
        if first[0] != '(':
            node['name'] = first[0:first.index('(')]
            pass

    ## if this node has a dist
        if last[-1] != ')':
            node['dist'] = float(last[last.rindex('):')+2:])
            pass

        return
    

    decorate(node, first, last, splittree)
    if node['leaf']: return node
    
    ## remove outer parens
    splittree[0] = first[first.index('(')+1:]
    splittree[-1] = last[:last.rindex(')')]
    paren_counts[0] -= 1
    paren_counts[-1] += 1
    

    for comma in range(1,len(paren_counts)):
        if debug: print 'considering comma', comma
        #if sum(paren_counts[:comma]) == 0 and sum(paren_counts[comma:]) == 0:
        if sum(paren_counts[:comma]) == 0:
            node['left'] = read_tree2('', splittree = splittree[:comma], paren_counts = paren_counts[:comma])
            node['right'] = read_tree2('', splittree = splittree[comma:], paren_counts = paren_counts[comma:])
            break
            pass
        pass
    
    node['left']['parent'] = node
    node['right']['parent'] = node
    
    return node


def dist_from_str(d, quick = True):
    # print 'dist_from_str', d

    if quick:
        if d == '': return 0
        return float(d)
    
    try:
        i = int(d)
        # print '  is int', i
        return i
    except ValueError:
        pass
    try:
        f = float(d)
        # print '  is float', f
        return f
    except ValueError:
        pass
    # print '  is nothing', 0
    return 0


def read_tree3(tree):
    
    rnode = None
    pnode = None
    in_dist = False
    ready_for_subnode = True

    # print tree

    for c in tree:
        still_in_dist = False
        # print 'considering', c, 'in dist:', in_dist, 'still in dist:', still_in_dist, 'ready for subnode:', ready_for_subnode, rnode

        if c == ' ':
            pass

        elif c == '(':
            ready_for_subnode = True
            node = {'name' : '',
                    'leaf' : False,
                    'parent' : pnode,
                    'dist' : '',
                    'children' : []}
            if pnode != None:
                pnode['children'].append(node)
                pass
            pnode = node
            
        elif c == ',':
            ready_for_subnode = True
            # close out the lower level node
            node['dist'] = dist_from_str(node['dist'])
            #node['dist'] = float(node['dist']) if node['dist'] != '' else 0

        elif c == ')':
            # close out the lower level node
            node['dist'] = dist_from_str(node['dist'])
            # fix left and right for pnode
            pnode['left'] = pnode['children'][0]
            pnode['right'] = pnode['children'][1]

            # and pop back up a level
            node = pnode
            pnode = node['parent']

        elif c == ';':
            # do nothing?
            pass

        elif c == ':':
            in_dist = True
            still_in_dist = True

        elif in_dist:
            node['dist'] += c
            still_in_dist = True

        elif ready_for_subnode:
            node = {'name' : c,
                     'leaf' : True,
                     'parent' : pnode,
                     'dist' : ''}
            if pnode != None:
                pnode['children'].append(node)
                pass
            ready_for_subnode = False

        else:
            node['name'] += c
            pass

        # this is how to get out of recording the distance (an of the key tokens will trigger this)
        if not still_in_dist: in_dist = False

        # save the root
        if rnode == None:
            rnode = node
            pass

        # print 'considered ', c, 'in dist:', in_dist, 'still in dist:', still_in_dist, 'ready for subnode:', ready_for_subnode, rnode

        pass

    # close out the root node
    rnode['dist'] = dist_from_str(rnode['dist'])

    return node
    

def new_inner_node(pnode):
    node = {'name' : '',
            'leaf' : False,
            'parent' : pnode,
            'dist' : '',
            'children' : []}
    if pnode != None:
        pnode['children'].append(node)
        pass
    return node

def new_leaf_node(pnode, c):
    node = {'name' : c,
            'leaf' : True,
            'parent' : pnode,
            'dist' : ''}
    if pnode != None:
        pnode['children'].append(node)
        pass
    return node

def read_tree4(tree, debug = False, contains = None, max_height = None):
    
    rnode = None
    pnode = None
    in_dist = False
    ready_for_subnode = True

    # print tree

    #for c in re.split('([():;,])', tree):
    for c in tree:
        still_in_dist = False
        if debug: print 'considering', c, 'in dist:', in_dist, 'still in dist:', still_in_dist, 'ready for subnode:', ready_for_subnode, rnode
        
        if c == ' ':
            continue

        elif c == '(':
            ready_for_subnode = True
            node = new_inner_node(pnode)
            pnode = node
            
        elif c == ',':
            ready_for_subnode = True
            # close out the lower level node
            node['dist'] = dist_from_str(node['dist'])
            # check for 'contains' node
            # if node['leaf'] and contains != None and node['name'] == contains:
            #     node['subtree_contains_node'] = True
                

        elif c == ')':
            # close out the lower level node
            node['dist'] = dist_from_str(node['dist'])
            # fix left and right for pnode
            pnode['left'] = pnode['children'][0]
            pnode['right'] = pnode['children'][1]

            # and pop back up a level
            node = pnode
            pnode = node['parent']

        elif c == ';':
            # do nothing?
            pass

        elif c == ':':
            in_dist = True
            still_in_dist = True

        elif in_dist:
            node['dist'] += c
            still_in_dist = True

        elif ready_for_subnode:
            node = new_leaf_node(pnode, c)
            ready_for_subnode = False

        else:
            node['name'] += c
            pass

        # this is how to get out of recording the distance (an of the key tokens will trigger this)
        if not still_in_dist: in_dist = False

        # save the root
        if rnode == None:
            rnode = node
            pass

        # print 'considered ', c, 'in dist:', in_dist, 'still in dist:', still_in_dist, 'ready for subnode:', ready_for_subnode, rnode

        pass

    # close out the root node
    rnode['dist'] = dist_from_str(rnode['dist'])

    return node
    






def find_node(tree, name):

    #print 'visiting', tree['name'], tree['name'] == name

    if tree['name'] == name:
        return tree

    if tree['leaf']:
        return None

    l = find_node(tree['left'], name)
    if l != None:
        return l

    r = find_node(tree['right'], name)
    if r != None:
        return r

    return None


def tree_to_str(node):
    
    if not node['leaf']:
        ret = '(%s, %s)' % (tree_to_str(node['left']), tree_to_str(node['right']))
    else:
        ret = ''
        pass

    if 'name' in node:
        ret = '%s%s' % (ret, node['name'])
        pass

    # print 'printing dist', node['name'], node['dist'], float(node['dist']) == 0.0, isinstance(node['dist'], int)
    if 'dist' in node and node['dist'] == 0:
        pass
    elif 'dist' in node and isinstance(node['dist'], int):
        ret = '%s:%d' % (ret, node['dist'])
    elif 'dist' in node:
        ret = '%s:%f' % (ret, node['dist'])
        pass
    
    return ret


def get_terminals(node):
    
    if node['leaf']:
        return [node]

    return get_terminals(node['left']) + get_terminals(node['right'])


## returns nodes in decending order
def get_path_to_root(node):
    if node['parent'] == None:
        return [node]
    return get_path_to_root(node['parent']) + [node]

def get_dist_btwn_nodes(ancestor, descendant):

    if ancestor == descendant: return 0

    if descendant['parent'] == None:
        print 'error in get_dist_btwn_nodes; ancestor not an ancestor of descendant?'
        print tree_to_str(ancestor)
        print tree_to_str(descendant)
        sys.exit(-1)
        pass

    return get_dist_btwn_nodes(ancestor, descendant['parent']) + descendant['dist']


if __name__ == "__main__":
    
    
    print '----------------------------'
    s = 'A:1'
    t = read_tree3(s)
    t2 = read_tree4(s)
    print s
    print t
    print t2
    print tree_to_str(t)
    print tree_to_str(t2)
    print find_node(t, 'B')
    print
    print '----------------------------'
    s = '(A, B):1;'
    t = read_tree3(s)
    t2 = read_tree4(s)
    print s
    print t
    print t2
    print tree_to_str(t)
    print tree_to_str(t2)
    print find_node(t, 'B')
    print
    print '----------------------------'
    s = '(A:1, (B:2, C:3)x:4)y:5;'
    t = read_tree3(s)
    t2 = read_tree4(s)
    print s
    print t
    print t2
    print tree_to_str(t)
    print tree_to_str(t2)
    print find_node(t, 'B')
    print
    print '----------------------------'
    s = '((D:2, A:1)yup, (B:2, C:3.0):4.1):5;'
    t = read_tree3(s)
    t2 = read_tree4(s)
    print s
    print t
    print t2
    print tree_to_str(t)
    print tree_to_str(t2)
    print find_node(t, 'B')
    

# sys.path.append('/net/akey/vol1/home/bvernot/archaic_exome/experiments/fdr_simulated_basic/latest/newick-1.3/build/lib/newick/')

# import newick
# print newick.parse_tree("(A,B);")



# class BranchLengthSum(newick.AbstractHandler):
#     def __init__(self):
#         self.sum = 0.0

#     def new_edge(self,b,l):
#         if l:
#             self.sum += l
#             pass

#     def get_result(self):
#         return self.sum

# print newick.parse('(A:1, (B:2, C:3):4):5;', BranchLengthSum())

