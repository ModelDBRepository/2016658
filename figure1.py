#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import neuron
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

h = neuron.h

dt = 0.01 #25

soma = h.Section()
soma.insert('pas')
soma.insert('nas')
soma.insert('kv3')
soma.insert('kv1')

# soma parameters
soma.diam = 20
soma.L = 126 #226
soma.Ra = 100.0
soma.cm = 1.0 # membrane capacitance (muF/cm2)

surface = soma.L * (soma.diam * 0.5)**2 * np.pi

# channel parameters
# passive resistance and resting Vm
soma(0.5).g_pas = 0.00025
soma(0.5).e_pas = -65
soma(0.5).gna_nas=0.1125
soma(0.5).thetam_nas=-22


# a-current
soma(0.5).gbar_kv1 = 0.005 
a_current_tau_scale = 7.5
neuron.h('a0h_kv1 = ' + str(a_current_tau_scale))
soma(0.5).ek = -90 # (mV)

h.psection(sec=soma)
print('A-current tau scale = ' + str(a_current_tau_scale))

# square pulse with 'IPSP' ramp
delay = 100.0
#max_amplitude = 0.01 #0.00055 # 0.02
reduction = 0 # 0.5 #0.5 # between 0 and 1
recovery = 200 # ramp duration
pre_duration = 300.0 #1000.0
post_duration = 00.0
min_amplitude = 0
backbaseline_duration = 100

fig = plt.figure(1,figsize=(5,12))
#fig size 5,12 to compare in fig 1
current=[]
rates=[]
#for fi curves use range(0,41) and stim_amp*0.025
#for traces for figure 1 use 0,13 and 130+0.03*stim_amp
for stim_amp in range(0,34): #35
    max_amplitude=stim_amp*0.03

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
    stim_electrode2 = h.IClamp(soma(0.5))
    stim_electrode2.dur = 1e9
    t_vec_stim = h.Vector([i * dt for i in range(len(stim_vec))])
    stim_vec.play(stim_electrode2._ref_amp, t_vec_stim, 1)


    # current_density = stim_electrode.amp / 1e3 / (surface / 1e8)

    I_record = h.Vector()
    Vm_record = h.Vector()
    Ia_record = h.Vector()


    I_record.record(stim_electrode2._ref_i)
    Vm_record.record(soma(0.5)._ref_v) #_ref_h_kv1 _ref_v
    Ia_record.record(soma(0.5)._ref_ik_kv1)

    tVec = h.Vector()
    tVec.record(h._ref_t)

    neuron.h.load_file('stdrun.hoc')
    neuron.h.dt = dt

    # Temperature (NEURON default = 6.3)
    temperature = 24
    neuron.h.celsius = temperature
    print('temperature = ' + str(temperature) + ' C')

    # initial Vm
    v_init = -65 # mV
    neuron.h('v_init=' + str(v_init))
    print('resting membrane potential = ' + str(v_init) + ' mV')
    neuron.h('init()')

    # Duration of simulation (in ms)
    neuron.run(1100.0)
    
    ax1 = fig.add_subplot(34, 1, stim_amp+1)
    ax1.plot(tVec, Vm_record)
    pA=round(max_amplitude*1000)
    ax1.set_ylabel(str(pA) +' pA')
    ax1.set_xlim(0,800)
    ax1.set_ylim(-100,60)
    #ax1.set_ylim(0.075,0.1)
    
    trace=Vm_record
    ti=tVec
    spikes=[]
    for x in range(0,(len(trace)-1),1):
        if trace[x]>=-30:
            if trace[x-1]<-30:
                        spikes.append(ti[x])
    if len(spikes)>0:
        if spikes[-1]>400:
                last_10=spikes[-10:-1]
                    #remove stuttering cells
                if len(np.diff(spikes))>0:
                        if max(np.diff(spikes)) <= 1.5*np.mean(np.diff(spikes)):
                            rates.append(1/np.mean(0.001*np.diff(last_10)))
                            current.append(max_amplitude)



plt.savefig("sim_to_traces.pdf",transparent=True)
#the trace to compare to is in the 20181001 folder abd file 18o01024
plt.show()
# print current_density


# In[2]:


plt.figure(figsize=(5, 5))
plt.plot(current, rates, color='C0',linestyle="-",marker="o")
plt.ylim(0,150)
plt.xlim(0,1)
plt.savefig("ficurve.png")
plt.show() 


# In[ ]:




