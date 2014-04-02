#https://gist.github.com/istepanov/6506508


"""
Search longest common substrings using generalized suffix trees built with Ukkonen's algorithm
 
Author: Ilya Stepanov <code at ilyastepanov.com>
 
(c) 2013
"""
 
 
#from __future__ import print_function
from Bio import SeqIO
from StringIO import StringIO
import sys
import re
import argparse


#Building Suffix Tree
 
END_OF_STRING = sys.maxint
 
 
class SuffixTreeNode:
    """
    Suffix tree node class. Actually, it also respresents a tree edge that points to this node.
    """
    new_identifier = 0
 
    def __init__(self, start=0, end=END_OF_STRING):
        self.identifier = SuffixTreeNode.new_identifier
        SuffixTreeNode.new_identifier += 1
 
        # suffix link is required by Ukkonen's algorithm
        self.suffix_link = None
 
        # child edges/nodes, each dict key represents the first letter of an edge
        self.edges = {}
 
        # stores reference to parent
        self.parent = None
 
        # bit vector shows to which strings this node belongs
        self.bit_vector = 0
 
        # edge info: start index and end index
        self.start = start
        self.end = end
 
    def add_child(self, key, start, end):
        """
        Create a new child node
 
        Agrs:
            key: a char that will be used during active edge searching
            start, end: node's edge start and end indices
 
        Returns:
            created child node
 
        """
        child = SuffixTreeNode(start=start, end=end)
        child.parent = self
        self.edges[key] = child
        return child
 
    def add_exisiting_node_as_child(self, key, node):
        """
        Add an existing node as a child
 
        Args:
            key: a char that will be used during active edge searching
            node: a node that will be added as a child
        """
        node.parent = self
        self.edges[key] = node
 
    def get_edge_length(self, current_index):
        """
        Get length of an edge that points to this node
 
        Args:
            current_index: index of current processing symbol (usefull for leaf nodes that have "infinity" end index)
        """
        return min(self.end, current_index + 1) - self.start
 
    def __str__(self):
        return 'id=' + str(self.identifier)
 
 
class SuffixTree:
    """
    Generalized suffix tree
    """
 
    def __init__(self):
        # the root node
        self.root = SuffixTreeNode()
 
        # all strings are concatenaited together. Tree's nodes stores only indices
        self.input_string = ''
 
        # number of strings stored by this tree
        self.strings_count = 0
 
        # list of tree leaves
        self.leaves = []
 
    def append_string(self, input_string):
        """
        Add new string to the suffix tree
        """
        start_index = len(self.input_string)
        current_string_index = self.strings_count
 
        # each sting should have a unique ending
        input_string += '$' + str(current_string_index)
 
        # gathering 'em all together
        self.input_string += input_string
        self.strings_count += 1
 
        # these 3 variables represents current "active point"
        active_node = self.root
        active_edge = 0
        active_length = 0
 
        # shows how many
        remainder = 0
 
        # new leaves appended to tree
        new_leaves = []
 
        # main circle
        for index in range(start_index, len(self.input_string)):
            previous_node = None
            remainder += 1
            while remainder > 0:
                if active_length == 0:
                    active_edge = index
 
                if self.input_string[active_edge] not in active_node.edges:
                    # no edge starting with current char, so creating a new leaf node
                    leaf_node = active_node.add_child(self.input_string[active_edge], index, END_OF_STRING)
 
                    # a leaf node will always be leaf node belonging to only one string
                    # (because each string has different termination)
                    leaf_node.bit_vector = 1 << current_string_index
                    new_leaves.append(leaf_node)
 
                    # doing suffix link magic
                    if previous_node is not None:
                        previous_node.suffix_link = active_node
                    previous_node = active_node
                else:
                    # ok, we've got an active edge
                    next_node = active_node.edges[self.input_string[active_edge]]
 
                    # walking down through edges (if active_length is bigger than edge length)
                    next_edge_length = next_node.get_edge_length(index)
                    if active_length >= next_node.get_edge_length(index):
                        active_edge += next_edge_length
                        active_length -= next_edge_length
                        active_node = next_node
                        continue
 
                    # current edge already contains the suffix we need to insert.
                    # Increase the active_length and go forward
                    if self.input_string[next_node.start + active_length] == self.input_string[index]:
                        active_length += 1
                        if previous_node is not None:
                            previous_node.suffix_link = active_node
                        previous_node = active_node
                        break
 
                    # splitting edge
                    split_node = active_node.add_child(
                        self.input_string[active_edge],
                        next_node.start,
                        next_node.start + active_length
                    )
                    next_node.start += active_length
                    split_node.add_exisiting_node_as_child(self.input_string[next_node.start], next_node)
                    leaf_node = split_node.add_child(self.input_string[index], index, END_OF_STRING)
                    leaf_node.bit_vector = 1 << current_string_index
                    new_leaves.append(leaf_node)
 
                    # suffix link magic again
                    if previous_node is not None:
                        previous_node.suffix_link = split_node
                    previous_node = split_node
 
                remainder -= 1
 
                # follow suffix link (if exists) or go to root
                if active_node == self.root and active_length > 0:
                    active_length -= 1
                    active_edge = index - remainder + 1
                else:
                    active_node = active_node.suffix_link if active_node.suffix_link is not None else self.root
 
        # update leaves ends from "infinity" to actual string end
        for leaf in new_leaves:
            leaf.end = len(self.input_string)
        self.leaves.extend(new_leaves)
 
    def find_longest_common_substrings(self):
        """
        Search longest common substrings in the tree by locating lowest common ancestors what belong to all strings
        """
 
        # all bits are set
        success_bit_vector = 2 ** self.strings_count - 1
 
        lowest_common_ancestors = []
 
        # going up to the root
        for leaf in self.leaves:
            node = leaf
            while node.parent is not None:
                if node.bit_vector != success_bit_vector:
                    # updating parent's bit vector
                    node.parent.bit_vector |= node.bit_vector
                    node = node.parent
                else:
                    # hey, we've found a lowest common ancestor!
                    lowest_common_ancestors.append(node)
                    break
 
        longest_common_substrings = ['']
        longest_length = 0
 
        # need to filter the result array and get the longest common strings
        for common_ancestor in lowest_common_ancestors:
            common_substring = ''
            node = common_ancestor
            while node.parent is not None:
                label = self.input_string[node.start:node.end]
                common_substring = label + common_substring
                node = node.parent
            # remove unique endings ($<number>), we don't need them anymore
            common_substring = re.sub(r'(.*?)\$?\d*$', r'\1', common_substring)
            if len(common_substring) > longest_length:
                longest_length = len(common_substring)
                longest_common_substrings = [common_substring]
            elif len(common_substring) == longest_length and common_substring not in longest_common_substrings:
                longest_common_substrings.append(common_substring)
 
        return longest_common_substrings
 
    def to_graphviz(self, node=None, output=''):
        """
        Show the tree as graphviz string. For debugging purposes only
        """
        if node is None:
            node = self.root
            output = 'digraph G {edge [arrowsize=0.4,fontsize=10];'
 
        output +=\
            str(node.identifier) + '[label="' +\
            str(node.identifier) + '\\n' + '{0:b}'.format(node.bit_vector).zfill(self.strings_count) + '"'
        if node.bit_vector == 2 ** self.strings_count - 1:
            output += ',style="filled",fillcolor="red"'
        output += '];'
        if node.suffix_link is not None:
            output += str(node.identifier) + '->' + str(node.suffix_link.identifier) + '[style="dashed"];'
 
        for child in node.edges.values():
            label = self.input_string[child.start:child.end]
            output += str(node.identifier) + '->' + str(child.identifier) + '[label="' + label + '"];'
            output = self.to_graphviz(child, output)
 
        if node == self.root:
            output += '}'
 
        return output
 
    def __str__(self):
        return self.to_graphviz()
 
# Read fasta files
def read_sequence(file_name):
    """
    Read the fasta file and return a string
    """
    handle = open(file_name)
    
    for record in SeqIO.parse(handle, "fasta"):
        #print record.id
        output = str(record.seq)
    handle.close()

        
    return output

# Recurively get all common substrings using suffix tree
def get_common_substrings(str_a, str_b, result, counter):

    if counter != 0:
        input_data = [str_a, str_b]
        suffix_tree = SuffixTree()
        for i in input_data:
            suffix_tree.append_string(i)

        lcs = suffix_tree.find_longest_common_substrings()

        for i in lcs:
            result.append(i)
            #print i
        
        for i in lcs:
            new_string_a = str_a.replace(i, "*")
            new_string_b = str_b.replace(i, "#")
            counter = counter - 1
        
        print "Counter is ", counter
        #print result
        get_common_substrings(new_string_a, new_string_b, result, counter)

    else:
        #print result
        return result


    
    
    
    
 
if __name__ == "__main__":

    string1 = read_sequence("sequenceA.fasta")
    string2 = read_sequence("sequenceB.fasta")
    print "Sequence A is....", string1
    print "Sequence B is....", string2
    print "Now trying to get all common substrings"
    anchor_candidates = []


    counter = 3
    get_common_substrings(string1, string2, anchor_candidates, counter)

    print anchor_candidates
##    for i in get_common_substrings(string1, string2, anchor_candidates):
##        print i


    
    
##    wait = input("PRESS ENTER TO CONTINUE")
##    suffix_tree = SuffixTree()
##    test_data = ["GATTTTAKBCCCC", "GATTTTZKBCCCC"]
##    test_data = [string1, string2]
##    for s in test_data:
##        suffix_tree.append_string(s)
##    
##    lcs = suffix_tree.find_longest_common_substrings()
##    print lcs
##    print "Printing..."
##    for s in lcs:
##        print s
##
##    new_string1 = string1.replace(lcs[0], "*")
##    new_string2 = string2.replace(lcs[0], "*")
##
##    test_data = [new_string1, new_string2]
##    for s in test_data:
##        suffix_tree.append_string(s)
##
##    lcs = suffix_tree.find_longest_common_substrings()
##    print lcs

