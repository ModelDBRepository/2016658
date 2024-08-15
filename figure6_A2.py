#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import neuron
import matplotlib.pyplot as plt

h = neuron.h

#mechanims = ['kv1', 'hh_wbm']
#for mech in mechanims:
#    neuron.load_mechanisms(mech)

dt = 0.01 #25

soma = h.Section()
soma.insert('pas')
soma.insert('nas')
soma.insert('kv3')
soma.insert('kv1')

# soma parameters
soma.diam = 20
soma.L = 126
soma.Ra = 100.0
soma.cm = 1.0 # membrane capacitance (muF/cm2)

surface = soma.L * (soma.diam * 0.5)**2 * np.pi

# channel parameters
# passive resistance and resting Vm
soma(0.5).g_pas = 0.00025
soma(0.5).e_pas = -65
soma(0.5).thetam_nas=-22

#Mess with stuff!
#soma(0.5).gkv3_kv3 = 0

# a-current
soma(0.5).gbar_kv1 = 0.005 # 0.01 returns interruption
a_current_tau_scale = 7.5 
neuron.h('a0h_kv1 = ' + str(a_current_tau_scale))
soma(0.5).ek = -90 # (mV)

h.psection(sec=soma)
print('A-current tau scale = ' + str(a_current_tau_scale))

# square pulse with 'IPSP' ramp
delay = 100.0
max_amplitude = 0.45 #2 # 
reduction = 0.222 # between 0 and 1
recovery = 200 # ramp duration
pre_duration = 1000.0
post_duration = 25000.0 #2000.0
min_amplitude = 0
backbaseline_duration = 1000

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
stim_electrode2 = h.IClamp(soma(0.5))
stim_electrode2.dur = 1e9
t_vec_stim = h.Vector([i * dt for i in range(len(stim_vec))])
stim_vec.play(stim_electrode2._ref_amp, t_vec_stim, 1)

noise_mean, noise_std_dev = 0, 0.16474

tvec = h.Vector(np.linspace(0.0, 25300, int(25300/0.01)))
noise_current = np.random.normal(noise_mean, noise_std_dev, len(tvec))
noise_current_vector = h.Vector(noise_current)

stim_electrode3 = h.IClamp(soma(0.5))
stim_electrode3.dur = 25300
#t_vec_stim = h.Vector([i * dt for i in range(len(noise_current))])
noise_current_vector.play(stim_electrode3._ref_amp, tvec, 1)
#################################################

##########

# current_density = stim_electrode.amp / 1e3 / (surface / 1e8)

I_record = h.Vector()
Vm_record = h.Vector()
Ia_record = h.Vector()
Ina_record = h.Vector()
Mna_record = h.Vector()
Hna_record = h.Vector()
Nkv3_record = h.Vector()
Inoise_record = h.Vector()
qkv1_record = h.Vector()
Ikv3_record = h.Vector()

I_record.record(stim_electrode2._ref_i)
Vm_record.record(soma(0.5)._ref_v)
Ia_record.record(soma(0.5)._ref_ik_kv1)
Ina_record.record(soma(0.5)._ref_ina_nas)
Inoise_record.record(stim_electrode3._ref_i)
qkv1_record.record(soma(0.5)._ref_q_kv1)
Mna_record.record(soma(0.5)._ref_m_nas)
Hna_record.record(soma(0.5)._ref_h_nas)
Nkv3_record.record(soma(0.5)._ref_n_kv3)
tVec = h.Vector()
tVec.record(h._ref_t)

neuron.h.load_file('stdrun.hoc')
neuron.h.dt = dt

# Temperature (NEURON default = 6.3)
temperature = 32
neuron.h.celsius = temperature
print('temperature = ' + str(temperature) + ' C')

# initial Vm
v_init = -65 # mV
neuron.h('v_init=' + str(v_init))
print('resting membrane potential = ' + str(v_init) + ' mV')
neuron.h('init()')

# Duration of simulation (in ms)
neuron.run(15500.0) #3500

##############################
# New section
# copy everything in this section to your code
# after neuron.run(...)
##############################
import os.path

# adjust this filename to where you want to store the traces
# careful: please leave the letter r in front of the quotes, 
# it is necessary so it works on Windows
output_name = input('Ask filename:')
r"C:output_name.txt"

# here I just make a quick check if the file already exists 
# so we don't overwrite it by accident
# if you want to overwrite it delete the old one first
#if os.path.exists(output_name):
#	e = "File already exists!"
#	raise RuntimeError(e)
with open(output_name, 'w') as output_file:
    for i, t in enumerate(tVec):
        out_line = '%f\t%f\t%f\t%f\n' % (t, Vm_record[i], I_record[i],Inoise_record[i])
        output_file.write(out_line)
# here we save the file as a simple .txt file
# first column time steps, second column Vm in mV, third column I in nA

##############################
# End new section
##############################

fig = plt.figure(1,figsize=(6,12))
ax1 = fig.add_subplot(7, 1, 1)
ax1.plot(tVec, Vm_record)
ax1.set_ylabel('Vm at soma (mV)')
#ax1.set_xlim(1000,2000)
#ax1.set_ylim(-50,-35)
ax2 = fig.add_subplot(7, 1, 2)
ax2.plot(tVec, I_record)
ax2.set_ylabel('Injected current (nA)')
ax3 = fig.add_subplot(7, 1, 3)
ax3.plot(tVec, Hna_record)
ax3.set_ylabel('Hna')
ax4 = fig.add_subplot(7, 1, 4)
ax4.plot(tVec, Mna_record)
ax4.set_ylabel('M_na')
ax5 = fig.add_subplot(7, 1, 5)
ax5.plot(tVec, Inoise_record)
ax5.set_ylabel('Inoise')
ax6 = fig.add_subplot(7, 1, 6)
ax6.plot(tVec, qkv1_record)
ax6.set_ylabel('qkv1')
ax7 = fig.add_subplot(7, 1, 7)
ax7.plot(tVec, Nkv3_record)
ax7.set_ylabel('N-kv3')
ax7.set_xlabel('Time (ms)')

plt.savefig("one_example_model.png")
plt.show()
# print current_density


# In[2]:


#Do the above in a loop, to make sure they are representative
fig = plt.figure(1,figsize=(6,12))
for x in range(0,10):
  print(x)
  np.random.seed(x)
  noise_current = np.random.normal(noise_mean, noise_std_dev, len(tvec))
  print(noise_current[0:5,])
  noise_current_vector = h.Vector(noise_current)
  
  stim_electrode3 = h.IClamp(soma(0.5))
  stim_electrode3.dur = 25300
  noise_current_vector.play(stim_electrode3._ref_amp, tvec, 1)
    
  Vm_record = h.Vector()
  Inoise_record = h.Vector()
    
  Vm_record.record(soma(0.5)._ref_v)
  Inoise_record.record(stim_electrode3._ref_i)
    
  tVec = h.Vector()
  tVec.record(h._ref_t)
  
  neuron.h.load_file('stdrun.hoc')
  neuron.h.dt = dt

    # Temperature (NEURON default = 6.3)
  temperature = 32
  neuron.h.celsius = temperature
  print('temperature = ' + str(temperature) + ' C')

    # initial Vm
  v_init = -65 # mV
  neuron.h('v_init=' + str(v_init))
  print('resting membrane potential = ' + str(v_init) + ' mV')
  neuron.h('init()') 
  neuron.run(25500.0)

  ax1 = fig.add_subplot(10, 1, x+1)
  ax1.plot(tVec, Vm_record)
  #ax1.set_xlim(0,0.05)
    
  output_name="permut"+str(x)+".csv"
  with open(output_name, 'w') as output_file:
      for i, t in enumerate(tVec):
          out_line = '%f\t%f\t%f\t%f\n' % (t, Vm_record[i], I_record[i],Inoise_record[i])
          output_file.write(out_line)
  #ax1.set_ylabel('Vm at soma (mV)')

plt.savefig("model_runs.png")
plt.show()


# In[ ]:




