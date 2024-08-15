#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Lets load some packages
from netpyne import sim, specs
#!pip install matplotlib==3.4.2 --user
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')
from neuron import gui #this line is absolutely necessary if importing cells from hoc templates
from neuron import h
#from neuron.units import um, mV, ms, cm
import math as math #for the sin function
import numpy as np
from numpy import random, exp, sqrt, log, pi
from random import shuffle, sample
h.load_file("stdrun.hoc")
#Lets load in Guillem's parameters for the neurons
N=1 #Params for 100 Neurons
params = np.loadtxt("neurons/params.dat"); # Parameters: EL, Rin and tauv, gNaf, rat (=gKf/gNaf), thm1, thh2, thn1, and gkv7f, tha1
ELs, Rins, tauvs, gNafs, rats, thm1s, thh2s, thn1s, gkv7fs, tha1s = params[:,:N];
Caps = tauvs/Rins
#To get a length from net capacitence
#most literature assumes neurons have 1uf/cm^2
#this capacitence is in nF
SA=Caps*1e-3
#neuron does length in um so lets make it um^2
SAum=SA*100000000
#we're going to assume a diameter of 20um to get a radius of 10um
#SA of a cyclinder(not including ends)=2*pi*r*L
lengths=SAum/(20*pi)
# Active conductances
# Heterogeneous peak conductances for active currents. Above are loaded normalized to Golomb values times capacitance
#this gave me nanosiemens I need siemens/cm^2 for it's conductance, so I need to multiple by e-9 to get nanosiemens then
#divide by the surface area in cm^2
gNas=np.array(gNafs)*112.5e+3*Caps*1e-9/SA
gKv3s=np.array(gNafs)*np.array(rats)*225.e+3*Caps*1e-9/SA
gkv7s=np.array(gkv7fs)*225.e+3*Caps*1e-9/SA
#Leak conductance is inverse of input resistance
gLs=1/Rins
#But that's in Mega Ohms and I need siemens
gLs=gLs*1e-6
#And I need it per cm^2
gLs=gLs/SA
# Network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters
netParams.stimSourceParams['iclamp'] = {'type': 'IClamp', 'amp': 0.0, 'dur': 2.0e3, 'delay': 0}
#I'm not happy with how I'm doing this, I don't think its the best way
#Normally I'd be able to make a cell type in cellParams and then make a population of that cell type with popParams
#But that's difficult to do in netpyne with trying to make them have different conductances and such
#For now it's making each individual cell a population
#there should be a way to fix this with cell lists, but the documentation on that's a little too limited for me to figure it out
for i in range(N):
    soma = {'geom': {}, 'mechs': {}}  # soma properties
    soma['geom'] = {'diam': 20, 'L': lengths[i], 'cm': 1}
    soma['mechs']['pas'] = {'g': gLs[i], 'e': ELs[i]}
    leak_conduct=gLs[i]
    soma['mechs']['naG'] = {'gbar': gNas[i], 'thm1': thm1s[i], 'thh2':thh2s[i] }
    soma['mechs']['kv7'] = {'gbar': gkv7s[i], 'tha1': tha1s[i] }
    soma['mechs']['kv3'] = {'gbar': gKv3s[i], 'thn1': thn1s[i] }
    soma['threshold'] = -10.0
    FS_dict = {'secs': {'soma': soma}}
    netParams.cellParams['FS'+str(i)] = FS_dict  # add rule dict to list of cell property rules
    netParams.stimTargetParams['iclamp->FS'+str(i)] = {'source': 'iclamp', 'conds': {'pop': 'FS'+str(i)}, 'sec': 'soma', 'loc': 0.5}


# In[2]:


dt = 0.01
# square pulse with 'IPSP' ramp
delay = 100.0
max_amplitude = 0.41 #0.41 #0.38 
reduction = 0.243 # between 0 and 1
recovery = 200 # ramp duration
pre_duration = 1000.0
post_duration = 2000.0
min_amplitude = 0
backbaseline_duration = 0

print('stim delay = ' + str(delay) + ' ms')
print('stim max amplitude = ' + str(max_amplitude) + ' nA')
print('stim reduction = ' + str(reduction))
print('stim recovery ramp = ' + str(recovery) + ' ms')
print('stim pre_duration = ' + str(pre_duration) + ' ms')
print('stim post_duration = ' + str(post_duration) + ' ms')

stim_amplitude = []
baseline_bins = int(delay / dt + 0.5)
for i in range(baseline_bins):
    stim_amplitude.append(0.0)
pre_duration_bins = int(pre_duration / dt + 0.5)
for i in range(pre_duration_bins):
    stim_amplitude.append(max_amplitude)
ipsp_duration_bins = int(recovery / dt + 0.5)
for i in range(ipsp_duration_bins):
    rel_duration = 1.0 * i / ipsp_duration_bins
    tmp_amplitude = (1 - reduction) * max_amplitude + rel_duration * reduction * max_amplitude
    stim_amplitude.append(tmp_amplitude)
post_duration_bins = int(post_duration / dt + 0.5)
for i in range(post_duration_bins):
    stim_amplitude.append(max_amplitude)
backbaseline_duration_bins = int(backbaseline_duration / dt + 0.5)
for i in range(backbaseline_duration_bins):
    stim_amplitude.append(0.0)
    
stim_vec = h.Vector(stim_amplitude)
t_vec_stim = h.Vector([i * dt for i in range(len(stim_vec))])

noise_mean, noise_std_dev = 0, 0.1*(0.41)

tvec = h.Vector(np.linspace(0.0, 3300, int(3300/0.01)))
noise_current = np.random.normal(noise_mean, noise_std_dev, len(tvec))
noise_current_vector = h.Vector()
noise_current_vector.from_python(noise_current)


# In[3]:


for i in range(N):
    netParams.popParams['CurrentClamp_pop'] = {'cellModel': 'FS'+str(i), 'cellType': 'FS'+str(i),  'numCells': 1}
    netParams.stimSourceParams['iclamp'] = {'type': 'IClamp',  'dur': 1e9, 'del':0}
    netParams.stimTargetParams['iclamp->CurrentClamp_pop'] = {
        'source': 'iclamp',
        'sec': 'soma',
        'loc': 0.5,
        'conds': {'pop':'CurrentClamp_pop'}}
    netParams.stimSourceParams['iclamp2'] = {'type': 'IClamp',  'dur': 1e9, 'del':0}
    netParams.stimTargetParams['iclamp2->CurrentClamp_pop'] = {
        'source': 'iclamp2',
        'sec': 'soma',
        'loc': 0.5,
        'conds': {'pop':'CurrentClamp_pop'}}
    simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration
    simConfig.duration = 3300
    simConfig.dt = 0.01  
    simConfig.recordStep = 0.01
    simConfig.verbose = False
    simConfig.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'},
                              'a_soma':{'sec':'soma','loc':0.5,'var':'a','mech':'kv7'},
                              'm_soma':{'sec':'soma','loc':0.5,'var':'m','mech':'naG'},
                              'h_soma':{'sec':'soma','loc':0.5,'var':'h','mech':'naG'},
                              'n_soma':{'sec':'soma','loc':0.5,'var':'n','mech':'kv3'}}
    simConfig.analysis['plotTraces'] = {'include':['all'],'saveFig': "interrupt/cell"+str(i)+".png"}#,'timeRange':[1100,1700]
    sim.create(netParams=netParams, simConfig=simConfig)
    for cell in sim.net.cells:
        stim_vec.play(cell.stims[0]['hObj']._ref_amp, t_vec_stim, True)
        #noise_current_vector.play(cell.stims[1]['hObj']._ref_amp, tvec, True)
    sim.simulate()
    sim.analyze()

    


# In[4]:


import pandas as pd
n=sim.allSimData['a_soma']['cell_0']
v=sim.allSimData['V_soma']['cell_0']

df = pd.DataFrame({'n':n, 'v':v}) #, 'n2':n2})

df.to_csv('nandvno_noise.csv')

df2 = pd.DataFrame({'i':stim_amplitude,'noise':noise_current,'t':(i * dt for i in range(len(stim_vec)))})
df2.to_csv("current_nonoise.csv")

df3 = pd.DataFrame({'i':stim_amplitude,'t':(i * dt for i in range(len(stim_vec)))})
df3.to_csv("current.csv")
#df3 = pd.DataFrame({'i':noise_current_vector,'t':(i * dt for i in range(len(noise_current_vector)))})
#df3.to_csv("currentnoise.csv")


# In[ ]:




