# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 21:22:28 2024

@author: liushiqing
"""

# for running PTSNet with KY networks
# PTSNet need to be installed first


# Import packages
import os
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from ptsnet.simulation.sim import PTSNETSimulation
from ptsnet.utils.io import get_example_path, export_time_series
import wntr as wn
import scipy.io as io
import pickle

#%% setting of the simulation
default_settings = {
    "time_step" : 0.001, # Simulation time step in [s]
    "duration" : 0.15, # Simulation duration in [s]
    "period" : 0, # Simulation period for EPSs
    # "default_wave_speed" : 1000, # Wave speed value for all pipes in [m/s]
    
    # select one of the following wave speed files
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky2/ky2_ws.txt', # Text file with wave speed values in KY2 network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky10/ky10_ws.txt', # Text file with wave speed values in KY10 network
    "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky15/ky15_ws.txt', # Text file with wave speed values in KY15 network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky18/ky18ws.txt', # Text file with wave speed values in KY18 network
    
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky2/ky2_ws_up.txt', # Text file with wave speed values in KY2_UP and KY2_UPIL network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky10/ky10_ws_up.txt', # Text file with wave speed values in KY10_UP and KY10_UPIL network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky15/ky15_ws_up.txt', # Text file with wave speed values in KY10_UP and KY15_UPIL network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky18/ky18ws_up.txt', # Text file with wave speed values in KY15_UP and KY18_UPIL network
    

    "delimiter" : '\t', # Delimiter of text file with wave speed values
    "wave_speed_method" : 'user', # Wave speed adjustment method
    "save_results" : False, # Saves numerical results in HDF5 format
    "skip_compatibility_check" : False, # Dismisses compatibility check
    "show_progress" : False, # Shows progress (Warnings should be off)
    "profiler_on" : False, # Measures computational times of the simulation
    "warnings_on" : False # Warnings are displayed if True
}

#%% test for a simple simulation without burst. 

#resutl matrix for no burst simulation

# select one of the following .inp files
# inpfile = 'C:/Users/liush/Desktop/KY/ky2/ky2_m.inp'                       # ky2 network with minor modifications
# inpfile = 'C:/Users/liush/Desktop/KY/ky10/ky10_m.inp'                   # ky10 network with minor modifications
inpfile = 'C:/Users/liush/Desktop/KY/ky15/ky15_m.inp'                   # ky15 network with minor modifications
# inpfile = 'C:/Users/liush/Desktop/KY/ky18/ky18_m.inp'                   # ky18 network with minor modifications

# inpfile = 'C:/Users/liush/Desktop/KY/ky2/ky2_m_up.inp'                     # ky2 network with unified pipes
# inpfile = 'C:/Users/liush/Desktop/KY/ky10/ky10_m_up.inp'                   # ky10 network with unified pipes
# inpfile = 'C:/Users/liush/Desktop/KY/ky15/ky15_m_up.inp'                   # ky15 network with unified pipes
# inpfile = 'C:/Users/liush/Desktop/KY/ky18/ky18_m_up.inp'                   # ky18 network with unified pipes

# inpfile = 'C:/Users/liush/Desktop/KY/ky2/ky2_m_upil.inp'                     # ky2 network with unified pipes and increased pipe lengths
# inpfile = 'C:/Users/liush/Desktop/KY/ky10/ky10_m_upil.inp'                   # ky10 network with unified pipes and increased pipe lengths
# inpfile = 'C:/Users/liush/Desktop/KY/ky15/ky15_m_upil.inp'                   # ky15 network with unified pipes and increased pipe lengths
# inpfile = 'C:/Users/liush/Desktop/KY/ky18/ky18_m_upil.inp'                   # ky18 network with unified pipes and increased pipe lengths



t1 = time.time()
simtest_noburst = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings)
numnode=simtest_noburst.wn.num_nodes # number of nodes in the network
nodelabel=simtest_noburst.ss['node'].labels
pipelabel=simtest_noburst.ss['pipe'].labels
simtest_noburst.run()

ky_noburst_sam=[]
pipev=[]
pipefr=[]

sample_feq=200 #frequency of sampling from original simulation result
ts=np.arange(0,default_settings['duration']/default_settings['time_step'],1/sample_feq/default_settings['time_step']) # timestamps of the target frequency
ts=ts.astype(np.int32)

for monitornodename in nodelabel[0:numnode]: 
    a=simtest_noburst['node'].head[str(monitornodename)]-simtest_noburst.ss['node'].elevation[str(monitornodename)]

    asam=a[ts] # only record the time stamps based on the sample_feq
    ky_noburst_sam=np.append(ky_noburst_sam,asam)

ky_noburst_sam_mat=np.reshape(ky_noburst_sam,(numnode,len(ts))) # reshape to a matrix


for monitorpipename in pipelabel[0:len(pipelabel)]:
    pipev_t =simtest_noburst.ss['pipe'].velocity[str(monitorpipename)] # record the velocity in pipes
    pipefr_t=simtest_noburst.ss['pipe'].flowrate[str(monitorpipename)] # record the flow rate in pipes
    
    pipev=np.append(pipev,pipev_t)
    pipefr=np.append(pipefr,pipefr_t)


elapsed = time.time() - t1
print(elapsed)


#%% find burst coefficients for target flow rate or pressure drop




simn = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings)
nodelabel=simn.ss['node'].labels # list of node label

startnode=20 # index of the start of burst node
endnode=30  # index of the end of burst node
nodeindex=np.arange(startnode,endnode)   # list of burst node indices

k=0
t0 = time.time()
bc0_temp=0.0003347 # initial guess of the burst coefficients
bc0=np.append([0]*startnode,[bc0_temp]*(endnode-startnode))    # the list of initial burst coefficients
lf0=[0]*startnode           # leak flow
pd0=[0]*endnode            # pressure drop


# run simulation with initial burst coefficients
for burstnodename in nodelabel[nodeindex]:
    simn = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings)    
    itemp=np.where(nodelabel == burstnodename)
    i=int(itemp[0])
    print(nodelabel[i])
    
    simn.add_burst(str(burstnodename), burst_coeff = bc0[i], start_time = 0.100, end_time = 0.115)
    # since we only care the pressure drop at source at t=0.115s, we can set the duration of simulation to slightly > 0.115s (e.g. 0.15s) to speed up this process
    simn.run()
    b=simn['node'].leak_flow[str(burstnodename)] # leak flow at burst node

    lf0=np.append(lf0,b[115]) # leak flow list
    
    a=simn['node'].head[str(burstnodename)] # pressure head at source node
    pd0[i]=a[100]-a[115] # pressure drop at source node

elapsedall = time.time() - t0
print(elapsedall)   

#%
#% find bc for target pressure drop/flow rate by iterations
lf=np.column_stack((lf0,lf0)) # duplicate leak flow 
bc=np.column_stack((bc0,bc0)) # duplicate burst coefficients
pdr=np.column_stack((pd0,pd0)) # duplicate pressure drop
pr=np.zeros(endnode) # pressure remaining

for ite in range(1,5): # run burst coefficients adjustment for 4 iterations
    print(ite)
    simn = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings)
    nodeindex=np.arange(startnode,endnode)   
    t0 = time.time()

    bc_temp1=5/pdr[startnode:endnode,ite]*bc[startnode:endnode,ite] # target pressure drop is 5m, estimate burst coefficient of the next iteration by linear calculation
    # bc_temp1=0.002/lf[:,ite]*bc[:,ite] # target: leak flow is 2L/s (0.002 m3/s), estimate burst coefficient of the next iteration by linear calculation
    
    bc_temp=np.append([0]*startnode,bc_temp1)    # the list of updated burst coefficients

    bc=np.column_stack((bc, bc_temp)) # update burst coefficients
    
    lftemp=[0]*startnode            # create a leak flow list
    pdtemp=[0]*startnode            # create a pressure drop list
    prtemp=[0]*startnode            # create a pressure remaining list
    
    for burstnodename in nodelabel[nodeindex]:
        simn = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings)    
        itemp=np.where(nodelabel == burstnodename)
        i=int(itemp[0])
        print(nodelabel[i])
        
        t1 = time.time()
        
        simn.add_burst(str(burstnodename), burst_coeff = bc_temp[i], start_time = 0.100, end_time = 0.115) # use the burst coefficient of the ith burst node
        simn.run()
        b=simn['node'].leak_flow[str(burstnodename)] # leak flow at current source node
        lftemp=np.append(lftemp,b[115]) # update leak flow vector

        a=simn['node'].head[str(burstnodename)] # pressure head at source node
        pdtemp=np.append(pdtemp,a[100]-a[115]) # pressure drop at source node
        prtemp=np.append(prtemp,a[115]-simn.ss['node'].elevation[i]) # pressure remaining
                    
    elapsedall = time.time() - t0
    print(elapsedall)    
    lf=np.column_stack((lf, lftemp)) # update leak flow matrix
    pdr=np.vstack((pdr.T,pdtemp)).T # update pressure drop matrix
    pr=np.vstack((pr.T,prtemp)).T # update pressure remaining matrix
    
#%% update the simulation duration
default_settings = {
    "time_step" : 0.001, # Simulation time step in [s]
    "duration" : 20, # Simulation duration in [s]
    "period" : 0, # Simulation period for EPSs
    # "default_wave_speed" : 1000, # Wave speed value for all pipes in [m/s]
    
    # select one of the following wave speed files
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky2/ky2_ws.txt', # Text file with wave speed values in KY2 network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky10/ky10_ws.txt', # Text file with wave speed values in KY10 network
    "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky15/ky15_ws.txt', # Text file with wave speed values in KY15 network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky18/ky18ws.txt', # Text file with wave speed values in KY18 network
    
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky2/ky2_ws_up.txt', # Text file with wave speed values in KY2_UP and KY2_UPIL network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky10/ky10_ws_up.txt', # Text file with wave speed values in KY10_UP and KY10_UPIL network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky15/ky15_ws_up.txt', # Text file with wave speed values in KY10_UP and KY15_UPIL network
    # "wave_speed_file_path" : 'C:/Users/liush/Desktop/KY/ky18/ky18ws_up.txt', # Text file with wave speed values in KY15_UP and KY18_UPIL network
    

    "delimiter" : '\t', # Delimiter of text file with wave speed values
    "wave_speed_method" : 'user', # Wave speed adjustment method
    "save_results" : False, # Saves numerical results in HDF5 format
    "skip_compatibility_check" : False, # Dismisses compatibility check
    "show_progress" : False, # Shows progress (Warnings should be off)
    "profiler_on" : False, # Measures computational times of the simulation
    "warnings_on" : False # Warnings are displayed if True
}    
    
#%% burst simulation at multiple sources


t1 = time.time()
simtest_noburst = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings) 
numnode=simtest_noburst.wn.num_nodes # number of nodes at which pressre are recorded
nodelabel=simtest_noburst.ss['node'].labels
simtest_noburst.run() # run burst free scenario

ky_noburst_sam=[]

sample_feq=200 #frequency of sampling from original simulation result
ts=np.arange(0,default_settings['duration']/default_settings['time_step'],1/sample_feq/default_settings['time_step'])
ts=ts.astype(np.int32) # time stamps which pressure is recorded

for monitornodename in nodelabel[0:numnode]: # record pressure
    a=simtest_noburst['node'].head[str(monitornodename)]-simtest_noburst.ss['node'].elevation[str(monitornodename)] # pressure at all time stamps at current monitoring node

    asam=a[ts] # pressure at the recorded time stampes
    ky_noburst_sam=np.append(ky_noburst_sam,asam) # append pressure vector of current monitoring node

ky_noburst_sam_mat=np.reshape(ky_noburst_sam,(numnode,len(ts))) # reshape to pressure matrix at burst-free scenario

simn = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings)

# Note: because the size of result file, only run no more than 150 simulations at a time
startnode=20
endnode=30  

nodeindex=np.arange(startnode,endnode)   
k=0
t0 = time.time()
bcfinal=bc[:,-1]    # use the burst coefficients at the final iteration
lf=[]            # leak flow
result_ky=list(range(endnode)) # result


for burstnodename in nodelabel[nodeindex]: # for each burst node
    simn = PTSNETSimulation(workspace_name = 'ky1', inpfile = inpfile, settings=default_settings)    
    itemp=np.where(nodelabel == burstnodename)
    i=int(itemp[0])
    print(nodelabel[i])
    
    t1 = time.time()
    
    simn.add_burst(str(burstnodename), burst_coeff = bcfinal[i], start_time = 0.100, end_time = 0.115)
    simn.run()
    b=simn['node'].leak_flow[str(burstnodename)] # record leak flow rate

    lf=np.append(lf,b[115]) # update leak flow vector
    
    result_sam=[] 

    for monitornodename in nodelabel[0:numnode]: # pressure monitoring node
        a=simn['node'].head[str(monitornodename)]-simn.ss['node'].elevation[str(monitornodename)] # pressure at all time stamps at current monitoring node
        asam=a[ts] # pressure at the recorded time stampes
        result_sam=np.append(result_sam,asam) # append pressure vector of current monitoring node
        
    result_sam=np.reshape(result_sam,(numnode,len(ts))) # reshape to pressure matrix

    r_mat_burstdiff=result_sam-ky_noburst_sam_mat[:,:] 
    # calculate the difference between the result of burst simulation and burst-free simulation, and the difference is the pressure change generated by induced bursts
    # This is to eliminate the effect of "numerical" transient, which is the pressure change occurred at the beginning of every simulation, including burst-free

    result_ky[nodeindex[k]]=r_mat_burstdiff # write into the result variable
    elapsed = time.time() - t1
    print(elapsed)
    k+=1
    
elapsedall = time.time() - t0
print(elapsedall)    

ky10_pd5m_0_100=result_ky # rename file
io.savemat("result_ky10_pd5m_0_100.mat",{'ky10_pd5m_0_100':ky10_pd5m_0_100}) # save file



    
    
    