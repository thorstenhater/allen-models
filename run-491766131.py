#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from collections import defaultdict
from pprint import pprint
import arbor
import json
import numpy as np

swc = 'Rbp4-Cre_KL100_Ai14-203503.04.01.01_527109145_m.swc'
fit = "491766131_fit.json"
currents = [0.150, 0.170, 0.210]

t_beg = 100
t_dur = 1000
t_end = t_beg + t_dur

def load_tree(swc, strip_axon=False, disable_taper=False):
    tree = arbor.load_swc(swc)
    if strip_axon:
        hd_axon, tl_axon = None, None
        ln_axon = 0
        x0, y0, z0 = None, None, None
        x1, y1, z1 = None, None, None
        new = arbor.sample_tree()
        for p, s in zip(tree.parents, tree.samples):
            if s.tag != 2:
                if (not tl_axon is None) and (p >= tl_axon):
                    p -= ln_axon
                new.append(parent=p, x=s.loc.x, y=s.loc.y, z=s.loc.z, radius=s.loc.radius, tag=s.tag)
            else:
                if hd_axon is None:
                    hd_axon = p
                    x0, y0, z0 = s.loc.x, s.loc.y, s.loc.z
                x1, y1, z1 = s.loc.x, s.loc.y, s.loc.z
                tl_axon = p
                ln_axon += 1
        axon = new.append(hd_axon, x0, y0, z0, 1, 2)
        new.append(axon, x0 + 60 + tree.samples[0].loc.radius, y0, z0, 1, 2)
        tree = new
    if disable_taper:
        new = arbor.sample_tree()
        for p, s in zip(tree.parents, tree.samples):
            new.append(parent=p, x=s.loc.x, y=s.loc.y, z=s.loc.z, radius=s.loc.radius, tag=s.tag)
            
    return tree

tree = load_tree(swc)

# Load mechanism data
with open(fit) as fd:
    fit = json.load(fd)
    ## collect parameters in dict
    mechs = defaultdict(dict)

### Passive parameters
ra = float(fit['passive'][0]['ra'])

### Remaining parameters
for block in fit['genome']:
    mech   = block['mechanism'] or 'pas'
    region = block['section']
    name = block['name']
    if name.endswith('_' + mech):
        name = name[:-(len(mech) + 1)]
    mechs[(mech, region)][name] = float(block['value'])

# Label regions
labels = arbor.label_dict({'soma': '(tag 1)',
                           'axon': '(tag 2)',
                           'dend': '(tag 3)',
                           'apic': '(tag 4)',
                           'center': '(location 0 0)'})

properties = fit['conditions'][0]
T  = properties['celsius'] + 273.15
Vm = properties['v_init']
    
# Run simulation
ts = []
volts = {}
counts = {}
for i in currents:
    morph = arbor.morphology(tree, spherical_root=True)

    # Build cell and attach Clamp and Detector
    cell = arbor.cable_cell(tree, labels)
    #cell.compartments_per_branch(1)
    cell.place('center', arbor.iclamp(100, 1000, 0.15))
    cell.place('center', arbor.spike_detector(-40))
    #cell.set_ion('ca', int_con=5e-5, ext_con=2.0, method=arbor.mechanism('nernst/x=ca'))
    cell.set_properties(tempK=T, Vm=Vm, rL=ra)

    print(' * Global Parameters')
    print(f"  * T  =  {T}K = {T - 273.15}C")    
    print(f"  * Vm =  {Vm}mV")    
    
    # Set reversal potentials
    print(" * Reversal Potential")
    for kv in properties['erev']:
        region = kv['section']
        for k, v in kv.items():
            if k == 'section':
                continue
            ion = k[1:]
            print(f'| {k:<3} | {region:6} | {v:10} |')
            cell.paint(region, arbor.ion(ion, rev_pot=float(v)))

    log_m = []
    log_p = []
            
    # Setup mechanisms and parameters
    print('Setting up mechanisms')
    ## Now paint the cell using the dict
    for (mech, region), vs in mechs.items():
        try:
            if mech != 'pas':
                m = arbor.mechanism(mech, vs)
                for k, v in vs.items():
                    log_p.append(f'| {mech:10} | {region:5} | {k:8} | {float(v):>12.6g} |')                
                log_m.append(f'||')
                cell.paint(region, m)
            else:
                m = arbor.mechanism('pas', {'e': vs['e'], 'g': vs['g']})
                for k, v in vs.items():
                    if k in {'e', 'g'}:
                        log_p.append(f'| {mech:10} | {region:5} | {k:8} | {float(v):>12.6g} |')                
                cell.paint(region, m)
                cell.paint(region, cm=vs["cm"]/100, rL=vs["Ra"])
        except Exception as e:
            print("ERROR")
            print("When trying to set", mech, vs)
            print("  ->", e)
            exit()            
    print('\n'.join(log_p))

        
    # Run the simulation, collecting voltages
    print('Simulation', end=' ')
    model = arbor.single_cell_model(cell)
    model.probe('voltage', 'center', frequency=10000)
    model.run(tfinal=1200)
    print('DONE')
    for t in model.traces:
        ts = t.time[:]
        volts[i] = t.value[:]
        break
    spikes = np.array(model.spikes)

    print("Simulation spikes:", spikes)
    count = len(spikes)
    print('Counted spikes', count)
    counts[i] = count
    break

# load references
ref = np.loadtxt('data.txt')
print(ref[0].shape)
# Plot voltage traces
fg, ax = plt.subplots()
shift = 900

#ax.bar(x=[t_beg], height=[200], width=[t_end - t_beg], bottom=[-90], align='edge', color='0.8', label='Stimulus')
for k, v in volts.items():
    ax.plot(np.array(ts), np.array(v) + 14, ls='-', label=f"I={k}mA Spikes={counts[k]}")
ax.plot(ref[0], ref[1] + 14, ls='-', label=f"Reference")
ax.set_ylabel('U/mV')
ax.set_xlabel('t/ms')
#ax.set_xlim(left=0, right=1200)
ax.set_ylim(bottom=-80, top=30)
ax.legend()
plt.show()
