python calibration_lib.py Tests/testTree.sim.lb.tre
normalize_tree.py Calibrated.tre Calibrated.norm.tre
lsd -c -i Deviated.tre -v  
nw_distance -n -mr -sa Deviated.tre.result.date.newick | sort > theirs
nw_distance -n -mr -sa Tests/testTree.sim.lb.norm.tre | sort > truth
nw_distance -n -mr -sa Calibrated.norm.tre | sort > mine
echo "sqrt(`paste truth mine theirs | awk '{print ($2-$4)*($2-$4)}' | numlist -avg`)" | bc
echo "sqrt(`paste truth mine theirs | awk '{print ($2-$6)*($2-$6)}' | numlist -avg`)" | bc
