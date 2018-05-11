from calibration_lib import calibrate_tree
from sys import argv
from dendropy import Tree

treefile = argv[1]
outfile = argv[2]

myTree = Tree.get_from_path(treefile,'newick')

print("Read tree successfully")

print(calibrate_tree(myTree,verbose=True))

myTree.write_to_path(outfile,"newick")


