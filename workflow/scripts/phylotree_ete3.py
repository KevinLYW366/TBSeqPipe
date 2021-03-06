# Use ETE3 to work with the newick tree which is generated by RAxML
## Reroot on canettii and prune

from ete3 import Tree
import sys

tree_in = sys.argv[1]
tree_out = sys.argv[2]

# read newick trees
t = Tree(tree_in, format = 0)

# reroot on canettii
t.set_outgroup("canettii")

# prune the tree to remove canettii and REF
node_list = []

for node in t:
    if node.name not in ["canettii", "REF"]:
        node_list.append(node)

t.prune(node_list)
t.write(format = 0, outfile = tree_out)



