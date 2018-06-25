#! /usr/bin/env python

from calibration_lib import calibrate_with_sampling_time,calibrate_log_opt,calibrate_composite_opt
from dendropy import TreeList
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-s","--samplingTime",required=True,help="Sampling time")
parser.add_argument("-r","--rootAge",required=False,help="Root age")
parser.add_argument("-t","--timeTree",required=True,help="The output trees with branch lengths in time unit")
parser.add_argument("-c","--composite",required=False,action='store_true',help="Do composite optimization. Default: NO")

args = vars(parser.parse_args())

myTrees = TreeList.get_from_path(args["input"],'newick')
smpl_times = {}
rootAge = float(args["rootAge"]) if args["rootAge"] else None

with open(args["samplingTime"],"r") as fin:
    fin.readline()
    for line in fin:
        name,time = line.split()
        smpl_times[name] = float(time)

for tree in myTrees:
    if args["composite"]:
        s = calibrate_composite_opt(tree,smpl_times,root_age=rootAge)
    else:    
        s = calibrate_log_opt(tree,smpl_times,root_age=rootAge,brScale=False)

myTrees.write_to_path(args["timeTree"],"newick")
print("Clock rate: " + str(s))
