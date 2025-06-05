"""
calc_relaxation.py
26 March 2025

Script to calculate N-H bond correlations from the MD trajectory. 

Incorporated from the GitHub repository: https://github.com/achicks15/CorrFunction_NMRRelaxation.git with a few modifications. 

NOTE: In main_implementation, cross check the names of the topology and trajectory file. 

python3 calc_relaxation.py <pdbpath>
"""


import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import mdtraj as md
import sys
import os

path = sys.argv[1]

print(path)

pdb=path.split("/")[-1]
# ff=sys.argv[2]
# B0 = sys.argv[3]
# params = sys.argv[4]

def calc_NHVecs(traj_file, top_file, start_snap=0, end_snap=-1):
    """
    Uses mdtraj to load the trajectory and get the atomic indices and coordinates to calculate the correlation functions.
    For each trajectory, load the trajectory using mdtraj, get the atomic index for the the N-H atoms and calculate the vector between the two.
    Append the vector to the NHVecs list for all the trajectories. 
    NHVecs should return a list of shape: (# Trajectories, # Snapshots, # Residues w/N-H Vectors, 3)
    """
    traj = md.load(traj_file, top=top_file)
    top = traj.topology
    
    # Process the trajectory such that the protein is aligned in all frames (remove overall translation and rotation)
    traj.superpose(traj, frame=0)

    # AtomSelection Indices
    Nit = top.select('name N and not resname PRO') ## PRO residue do not have N-H vectors
    Hyd = top.select('name H and not resname PRO')
    NH_Pair = [[i,j] for i,j in zip(Nit,Hyd)]
    NH_Pair_Name = [[top.atom(i),top.atom(j)] for i,j in NH_Pair]
    NH_Res = ["{}-{}{}".format(str(i).split('-')[0],str(i).split('-')[1], str(j).split('-')[1]) for i,j in NH_Pair_Name]
    # print(NH_Res)

    ##Generate the N-H vectors in Laboratory Frame
    NHVecs_tmp = np.take(traj.xyz, Hyd, axis=1) - np.take(traj.xyz, Nit, axis=1)
    sh = list(NHVecs_tmp.shape)
    sh[2] = 1
    NHVecs_tmp = NHVecs_tmp / np.linalg.norm(NHVecs_tmp, axis=2).reshape(sh)
    
    return NHVecs_tmp[start_snap:end_snap]

def split_NHVecs(nhvecs, dt, tau):
    """
    This function will split the trajectory in chunks defined by tau. 
    nhvecs = array of N-H bond vectors,
    dt = timestep of the simulation
    tau = length of chunks
    """
    nFiles = len(nhvecs) ## number of trajectories
    nFramesPerChunk = int(tau/dt) ###tau/timestep 
    print("Number of frames per chunk: ", nFramesPerChunk)
    used_frames = np.zeros(nFiles,dtype=int)
    remainingFrames = np.zeros(nFiles,dtype=int)
    for i in range(nFiles):
        nFrames = nhvecs[i].shape[0]
        used_frames[i] = int(nFrames/nFramesPerChunk)*nFramesPerChunk
        remainingFrames[i] = nFrames % nFramesPerChunk
    
    nFramesTot=int(used_frames.sum())
    out = np.zeros((nFramesTot,nhvecs[0].shape[1],nhvecs[0].shape[2]), dtype=nhvecs[0].dtype)
    start = 0
    for i in range(nFiles):
        end = int(start+used_frames[i])
        endv = int(used_frames[i])
        out[start:end,...] = nhvecs[i][0:endv,...]
        start = end
        
    sh = out.shape
    vecs = out.reshape((int(nFramesTot/nFramesPerChunk), nFramesPerChunk, sh[-2], sh[-1]))
    
    return vecs

def calc_Ct(nhvecs):
    """
    Calculates the correlation function of the N-H bond vectors found in nhvecs. 
    Direct space calculation. This could be changed to Fourier space calculation for increased speed. 
    
    LICENSE INFO:
    
    MIT License

    Copyright (c) 2017 Po-chia Chen

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """
    sh = nhvecs.shape
    print(sh)
    nReplicates=sh[0] ; nDeltas=int(sh[1]/2) ; nResidues=sh[2]
    print(nReplicates, nResidues, nDeltas)
    Ct  = np.zeros( (nDeltas, nResidues), dtype=nhvecs.dtype )
    dCt = np.zeros( (nDeltas, nResidues), dtype=nhvecs.dtype )
    
    for delta in range(1,1+nDeltas):
        nVals=sh[1]-delta
        # = = Create < vi.v'i > with dimensions (nRep, nFr, nRes, 3) -> (nRep, nFr, nRes) -> ( nRep, nRes ), then average across replicates with SEM.
        tmp = -0.5 + 1.5 * np.square( np.einsum( 'ijkl,ijkl->ijk', nhvecs[:,:-delta,...] , nhvecs[:,delta:,...] ) )
        tmp  = np.einsum( 'ijk->ik', tmp ) / nVals
        Ct[delta-1]  = np.mean( tmp, axis=0 )

        # Modified denominator to incorporate only 1 simulation replicate and avoid invalid divisions.
        dCt[delta-1] = np.std( tmp, axis=0 ) /  max(np.sqrt(nReplicates), 1.0)
        # dCt[delta-1] = np.std( tmp, axis=0 ) / ( np.sqrt(nReplicates) - 1.0 )
    
    return Ct, dCt

def main_implementation(pdb):

    ## Global Variables for the calculation of the NH Vecs and the correlation functions

    FileLoc = "" ## Main Directory Location
    FTOPN = pdb+"_nowater.gro" ## Name of topology for the trajectory
    FMDN = pdb+"_nowater.xtc"  ## Name of the trajectory, should be centered and stripped of solute
    # FTOPN = "md_centered.gro"
    # FMDN = "md_centered.xtc"

    ## Parameters and Physical Constants for calculation of Relaxation Rates

    # H_gyro = 2*np.pi*42.57748*1e6     ## Gyromagnetic Ratio: Hydrogen ([rad]/[s][T]) 
    # N_gyro = -2*np.pi*4.317267*1e6     ## Gyromagnetic Ratio: Nitrogen ([rad]/[s][T]) 

    ## Need 5 Frequencies: ## J[0], J[wH], J[wN], J[wH-wN], J[wH+wN]
    # Larmor1H = H_gyro*B0              ## Larmor Frequency: Hydrogen ([rad]/[s])
    # Larmor15N = N_gyro*B0             ## Larmor Frequency: Nitrogen ([rad]/[s])
    # omDiff = Larmor1H - Larmor15N    ## Diff in Larmor Frequencies of Spin IS
    # omSum  = Larmor1H + Larmor15N    ## Sum of Larmor Frequencies of Spin IS

    # mu_0 = 4*np.pi*1e-7    ; ## Permeability of Free Space: ([H]/[m]) 
    # hbar = 1.0545718e-34  ; ## Reduced Plank's constant: [J] * [s] = [kg] * [m^2] * [s^-1] 

    # R_NH = 1.02e-10                     ## distance between N-H atoms in Angstroms
    # dSigmaN = -170e-6               ##  CSA of the S-spin atom

    # FDD = (1./10.)*np.power((mu_0*hbar*H_gyro*N_gyro)/(4*np.pi*np.power(R_NH,3)),2)
    #FCSA = 498637299.69233465
    # FCSA = (2.0/15.0)*(Larmor15N**2)*(dSigmaN**2)        ## CSA factor 

    # Make a new directory to store relaxation time and correlation data
    os.chdir(path+'/Analysis_equi/processed_trajs/')
    # os.chdir(path)
    # if not os.path.isfile("../NMR_analysis/relaxation_data"):
        # os.mkdir("../NMR_analysis/relaxation_data")

    ## Calculate the NHVecs; Can be adapted to loop over multiple trajectories using TRAJLIST_LOC
    NHVecs = []
    start=0; end=-1;  ## 
    NHV = calc_NHVecs(FMDN, FTOPN, start, end)
    NHVecs.append(NHV)

    dt = 200 ## timestep of simulations: (ps)
    tau_split = np.array(NHVecs).shape[1]*dt ## Number of snapshots to calculate the correlation function over.

    print("Shape of NH vector array: ", np.array(NHVecs).shape[1])
    print("Simulation time step: " , dt)
    print("Number of snapshots to calculate correlation function over: ", tau_split)

    ## Split the vecs based off the tau_split you want and the time step. 
    vecs_split = split_NHVecs(NHVecs, dt, tau_split) 

    ## Calculate the correlation functions and the standard deviation in the correlation function.
    ## Save the correlation functions in a dataframe and then to a csv file for later use.
    Ct, dCt = calc_Ct(vecs_split)
    
    ## Convert to dataframe with index set as timesteps in ns
    filepath="../NMR_analysis/relaxation_data/"
    # print current working directory
    print("Current working directory: ", os.getcwd())

    CtOutFname = filepath+'NH_Ct.csv'
    dCtOutFname = filepath+'NH_dCt.csv'
    CtDF = pd.DataFrame(Ct, index = np.arange(1, Ct.shape[0]+1)*dt/1000) 
    dCtDF = pd.DataFrame(dCt, index = np.arange(1, dCt.shape[0]+1)*dt/1000)
    CtDF.to_csv(CtOutFname)
    dCtDF.to_csv(dCtOutFname)

main_implementation(pdb)
