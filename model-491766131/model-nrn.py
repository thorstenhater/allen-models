#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

import time
from collections import defaultdict
import json

import arbor

from allensdk.model.biophys_sim.config import Config
from allensdk.model.biophysical.utils import create_utils

def run_arb(fit, swc, current, t_start, t_stop):
    tree = arbor.load_swc(swc)

    # Load mechanism data
    with open(fit) as fd:
        fit = json.load(fd)
    ## collect parameters in dict
    mechs = defaultdict(dict)

    ### Passive parameters
    ra = float(fit['passive'][0]['ra'])

    ### Remaining parameters
    for block in fit['genome']:
        mech = block['mechanism'] or 'pas'
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
                               'center': '(location 0 0.5)'})

    properties = fit['conditions'][0]
    T  = properties['celsius'] + 273.15
    Vm = properties['v_init']

    # Run simulation
    morph = arbor.morphology(tree, spherical_root=True)

    # Build cell and attach Clamp and Detector
    cell = arbor.cable_cell(morph, labels)
    cell.place('center', arbor.iclamp(t_start, t_stop - t_start, current))
    cell.place('center', arbor.spike_detector(-40))
    cell.compartments_length(20)
    
    # read json file and proceed to set parameters and mechanisms

    # set global values
    print('Setting global parameters')
    print(f"  * T  =  {T}K = {T - 273.15}C")
    print(f"  * Vm =  {Vm}mV")
    cell.set_properties(tempK=T, Vm=Vm, rL=ra)

    # Set reversal potentials
    print("Setting reversal potential for")
    for kv in properties['erev']:
        region = kv['section']
        for k, v in kv.items():
            if k == 'section':
                continue
            ion = k[1:]
            print(f'  * region {region:6} species {ion:5}: {v:10}')
            cell.paint(region, arbor.ion(ion, rev_pot=float(v)))

    cell.set_ion('ca', int_con=5e-5, ext_con=2.0, method=arbor.mechanism('default_nernst/x=ca'))

    # Setup mechanisms and parameters
    print('Setting up mechanisms')
    ## Now paint the cell using the dict
    for (mech, region), vs in mechs.items():
        print(f"  * {region:10} -> {mech:10}: {str(vs):>60}", end=' ')
        try:
            if mech != 'pas':
                m = arbor.mechanism(mech, vs)
                cell.paint(region, m)
            else:
                m = arbor.mechanism('default_pas', {'e': vs['e'], 'g': vs['g']})
                cell.paint(region, m)
                cell.paint(region, cm=vs["cm"]/100, rL=vs["Ra"])
            print("OK")
        except Exception as e:
            print("ERROR")
            print("When trying to set", mech, vs)
            print("  ->", e)
            exit()

    # Run the simulation, collecting voltages
    print('Simulation', end=' ')
    default = arbor.default_catalogue()
    catalogue = arbor.allen_catalogue()
    catalogue.insert(default, 'default_')
    
    model = arbor.single_cell_model(cell)
    model.properties.catalogue = catalogue
    model.probe('voltage', 'center', frequency=200000)
    model.run(tfinal=t_start + t_stop, dt=1000/200000)
    print('DONE')
    for t in model.traces:
        ts = t.time[:]
        vs = t.value[:]
        break
    spikes = np.array(model.spikes)
    count = len(spikes)
    print('Counted spikes', count)
    return np.array(ts), np.array(vs) + 14

def run_nrn(manifest, current, t_start, t_stop):
    description = Config().load(manifest)

    description.update_data({'axon_type': 'truncated'}, 'biophys')

    fix_sections = ['passive', 'axon_morph,', 'conditions', 'fitting']
    description.fix_unary_sections(fix_sections)

    utils = create_utils(description)
    hoc   = utils.h

    manifest = description.manifest
    morphology_path = manifest.get_path('MORPHOLOGY').encode('ascii', 'ignore').decode("utf-8")
    stimulus_path = description.manifest.get_path('stimulus_path')

    utils.generate_morphology(morphology_path)
    utils.load_cell_parameters()
    utils.read_stimulus(stimulus_path, sweep=35) # needed for sampling rates

    stim       = hoc.IClamp(hoc.soma[0](0.5))
    stim.amp   = current
    stim.delay = t_start
    stim.dur   = t_stop - t_start

    simulation_dt = 1.0e3/utils.simulation_sampling_rate

    hoc.dt = simulation_dt
    hoc.tstop = t_start + t_stop

    vec = utils.record_values()
    hoc.finitialize()
    hoc.run()

    recorded_data = utils.get_recorded_data(vec)

    vs = recorded_data['v']
    ts = recorded_data['t']

    return ts, vs

manifest = 'manifest.json'
swc      = 'Rbp4-Cre_KL100_Ai14-203503.04.01.01_527109145_m.swc'
fit      = "491766131_fit.json"
current  = 0.15
t_start  = 200
t_stop   = 1200

g, ax = plt.subplots()
ax.bar(x=[t_start], height=[200], width=[t_stop - t_start], bottom=[-100], align='edge', color='0.8', label='Stimulus')

t0 = time.perf_counter()
ts, vs = run_nrn(manifest, current, t_start, t_stop)
t1 = time.perf_counter()
print(f"t_nrn={t1 - t0}")

ax.plot(ts, vs*1000, label='Allen', ls='-')

t2 = time.perf_counter()
ts, vs = run_arb(fit, swc, current, t_start, t_stop)
t3 = time.perf_counter()
print(f"t_arb={t3 - t2}")

ax.plot(ts, vs, ls='-', label='Arbor')

ax.set_ylabel('U/mV')
ax.set_xlabel('t/ms')
ax.set_xlim(left=0, right=t_start + t_stop)
ax.set_ylim(bottom=-90, top=10)
ax.legend()
plt.show()
