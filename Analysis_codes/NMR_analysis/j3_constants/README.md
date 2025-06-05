Computes J3 coupling constants from the trajectory
Once ```process_trajectories``` is run, just run the following to automatically parse J3 constants, compute them from the trajectory, integrate data, compute relative errors and plot them for each residue.

```
./run_j3_main.sh <path> <offset> <force field>
```
```offset``` is used to make sure that the same residues are mapped in the simulation trajectory and the experimental dataset.
