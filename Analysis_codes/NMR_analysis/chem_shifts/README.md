Computes chemical shifts from the trajectory using MD2NMR module. 
Once ```process_trajectories``` is run, just run the following to automatically parse experimental shifts, compute shifts from the trajectory, integrate data, compute relative errors and plot them for each residue.

```
./run_shifts_main.sh <path> <offset> <force field>
```
```offset``` is used to make sure that the same residues are mapped in the simulation trajectory and the experimental dataset.
Mention ```force field``` as either 'amber' or 'charmm' in this code, to choose the appropriate MD2NMR version. 
