#!/usr/bin/env python
# coding: utf-8

# In[1]:


### import numpy as np
import neuron
import matplotlib.pyplot as plt
import numpy as np
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
a_current_tau_scale =  7.5
neuron.h('a0h_kv1 = ' + str(a_current_tau_scale))
soma(0.5).ek = -90 # (mV)

h.psection(sec=soma)
print('A-current tau scale = ' + str(a_current_tau_scale))

# square pulse with 'IPSP' ramp
delay = 100.0
max_amplitude = 0.45 #0.00055 # 0.02
reduction = 0.222 #0.5 # between 0 and 1
recovery = 200 # ramp duration
pre_duration = 1000.0 #1000.0
post_duration = 2500.0
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



# In[ ]:


subplotnum=1
for tex in [1400.0, 1600.00, 1800.0, 2000.00, 2200.0]:
    for cond in range(10,200,10) :
        ##########################################################
        # model of excitatory and inhibitory synaptic input
        ###########################################################
        # parameters
        # max. amplitude for E synaptic input; adjust until appropriate
        e_syn_amp = cond*0.0001       
        # max. amplitude for I synaptic input; adjust until appropriate
        #i_syn_amp = 10.000           
        # reversal potentials in mV
        e_rev = -65.0 #0.0
        #i_rev = -65.0
        # rise and decay time constants; here I added typical values for AMPA and GABA-A receptors
        e_tau_rise = 0.5
        e_tau_decay = 1
        #i_tau_rise = 1.0
        #i_tau_decay = 20.0 #20 is default
        # activation time of E synapse in ms
        i_activation_time = tex

        # now we set up the two synapses for NEURON; no changes required here
        e_syn = h.Exp2Syn(soma(0.5))
        e_syn.tau1 = e_tau_rise
        e_syn.tau2 = e_tau_decay
        e_syn.e = e_rev
        i_stim_vec = h.Vector([i_activation_time])
        vs_i = h.VecStim()
        vs_i.play(i_stim_vec)
        nc_i = h.NetCon(vs_i, e_syn)
        nc_i.delay = 0.0
        nc_i.weight[0] = e_syn_amp
        nc_i.threshold = 0.0
        ###########################################################


        # current_density = stim_electrode.amp / 1e3 / (surface / 1e8)

        I_record = h.Vector()
        Vm_record = h.Vector()
        Ia_record = h.Vector()
        Ina_record = h.Vector()
        Mna_record = h.Vector()
        Hna_record = h.Vector()
        Nkv3_record = h.Vector()
        qkv1_record = h.Vector()
        qtau_record = h.Vector()
        pkv1_record = h.Vector()
        Ikv3_record = h.Vector()

        I_record.record(stim_electrode2._ref_i)
        Vm_record.record(soma(0.5)._ref_v)
        Ia_record.record(soma(0.5)._ref_ik_kv1)
        Ina_record.record(soma(0.5)._ref_ina_nas)
        qkv1_record.record(soma(0.5)._ref_q_kv1)
        pkv1_record.record(soma(0.5)._ref_p_kv1)
        Mna_record.record(soma(0.5)._ref_m_nas)
        Hna_record.record(soma(0.5)._ref_h_nas)
        Nkv3_record.record(soma(0.5)._ref_n_kv3)
        qtau_record.record(soma(0.5)._ref_qtau_kv1)
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
        neuron.run(4500.0)

        ##############################
        # New section
        # copy everything in this section to your code
        # after neuron.run(...)
        ##############################
        import os.path

        # adjust this filename to where you want to store the traces
        # careful: please leave the letter r in front of the quotes, 
        # it is necessary so it works on Windows
        output_name = "syn_runs_inht"+str(round(tex))+"cond"+str(cond)
        output_name=output_name.replace(".","_")+".txt"

        # here I just make a quick check if the file already exists 
        # so we don't overwrite it by accident
        # if you want to overwrite it delete the old one first
        #if os.path.exists(output_name):
        #	e = "File already exists!"
        #	raise RuntimeError(e)

        # here we save the file as a simple .txt file
        # first column time steps, second column Vm in mV, third column I in nA
        with open(output_name, 'w') as output_file:
            for i, t in enumerate(tVec):
                out_line = '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (t, Vm_record[i], I_record[i], Ia_record[i], Ina_record[i],
                                                                 qkv1_record[i], pkv1_record[i],Mna_record[i],Hna_record[i],Nkv3_record[i]
                                                                )
                output_file.write(out_line)

        ##############################
        # End new section
        ##############################

        fig = plt.figure(1,figsize=(6,28))
        ax1 = fig.add_subplot(20, 1, subplotnum)
        ax1.plot(tVec, Vm_record)
        ax1.set_ylabel(str(cond))
        ax1.set_xlim(1000,3000)
        subplotnum+=1
        
    plot_name="syn_runs_inh"+str(round(tex))+".png"
    plt.savefig(plot_name.replace(".","_"))
    subplotnum=1
        #plt.show()
        # print current_density






