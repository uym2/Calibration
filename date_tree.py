#! /usr/bin/env python

from calibration_lib import calibrate_with_sampling_time,calibrate_log_opt
from dendropy import TreeList
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-s","--samplingTime",required=True,help="Sampling time")
parser.add_argument("-r","--rootAge",required=True,help="Root age")
parser.add_argument("-t","--timeTree",required=True,help="The output trees with branch lengths in time unit")
#parser.add_argument("-d","--nodeDate",required=False,help="Output all node dates")

args = vars(parser.parse_args())

myTrees = TreeList.get_from_path(args["input"],'newick')
smpl_times = {}
rootAge = float(args["rootAge"]) if args["rootAge"] else None

with open(args["samplingTime"],"r") as fin:
    fin.readline()
    for line in fin:
        name,time = line.split()
        smpl_times[name] = float(time)

#with open(args["nodeDate"],'w') as fout:
#i = 0
for tree in myTrees:
    #f = calibrate_with_sampling_time(tree,smpl_times)
    f = calibrate_log_opt(tree,smpl_times,root_age=rootAge)
    #fout.write("Tree " + str(i) + "\n")
    #dates = {}
    #for node in tree.postorder_node_iter():
        #label = node.taxon.label if node.is_leaf() else node.label
        #dates[label] = node.time
    #for label in sorted(dates.keys()):
    #    fout.write(label + " " + str(dates[label]) + "\n")    
    #i += 1

myTrees.write_to_path(args["timeTree"],"newick")
