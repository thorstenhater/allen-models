#!/usr/bin/env python3

import matplotlib.pyplot as plt

from collections import defaultdict
from pprint import pprint
import arbor
import json
import numpy as np

ts = []
volts = {}
counts = {}
for i in [0.130, 0.190, 0.270]:
    # Load tree and build morphology
    tree = arbor.load_swc('473561729/Gad2-IRES-Cre_Ai14_IVSCC_-172679.03.01.01_471076778_m.swc')



    # copy out the contents
    p_axon = None
    x0, y0, z0 = None, None, None
    x1, y1, z1 = None, None, None
    new = arbor.sample_tree()
    for p, s in zip(tree.parents, tree.samples):
        if s.tag != 2:
            new.append(parent=p, x=s.loc.x, y=s.loc.y, z=s.loc.z, radius=s.loc.radius, tag=s.tag)
        else:
            if p_axon is None:
                p_axon = p
                x0, y0, z0 = s.loc.x, s.loc.y, s.loc.z
            x1, y1, z1 = s.loc.x, s.loc.y, s.loc.z

    axon = new.append(p_axon, x0, y0, z0, 1, 2)
    new.append(axon, x0 + 60 + tree.samples[0].loc.radius, y0, z0, 1, 2)

    morph = arbor.morphology(tree, spherical_root=True)

    # Label regions
    labels = arbor.label_dict({'soma': '(tag 1)',
                               'axon': '(tag 2)',
                               'dend': '(tag 3)',
                               'apic': '(tag 4)',
                               'center': '(location 0 0.5)'})

    # Build cell and attach Clamp and Detector
    cell = arbor.cable_cell(tree, labels)

    cell.place('center', arbor.iclamp(200, 1200, i))
    cell.place('center', arbor.spike_detector(-40))

    # read json file and proceed to set parameters and mechanisms
    with open('473561729/473561729_fit.json') as fd:
        fit = json.load(fd)

    # set global values
    print('Setting global parameters')
    properties = fit['conditions'][0]
    T  = properties['celsius'] + 273.15
    print(f"  * T  =  {T}K = {T - 273.15}C")
    Vm = properties['v_init']
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
            print(f'  * species {ion:5}: {v:10}')
            cell.paint(region, arbor.ion(ion, rev_pot=float(v)))

    # Setup mechanisms and parameters

    ## collect parameters in dict
    mechs = defaultdict(dict)

    ### Passive parameters
    pas = fit['passive'][0]
    e_pas  = pas['e_pas']
    ra_pas = pas['ra']#/100

    print('Setting passive parameters')
    for v in pas['cm']:
        cm_pas = float(v['cm'])/100
        region = v['section']
        print(f"  * {region:10} -> rL={ra_pas}, Vm={e_pas} cm={cm_pas}")
        cell.paint(region, cm=cm_pas, rL=ra_pas, Vm=e_pas,)

    ### Remaining parameters
    for block in fit['genome']:
        mech   = block['mechanism'] or 'pas'
        region = block['section']
        name   = block['name'][:-(len(mech) + 1)]
        mechs[(mech, region)][name] = block['value']

    print('Setting up mechanisms')
    ## Now paint the cell using the dict
    for (mech, region), vs in mechs.items():
        print(f"  * {region:10} -> {mech:10}: {str(vs):60}", end=' ')
        try:
            m = arbor.mechanism(mech, vs)
            cell.paint(region, m)
            print("OK")
        except Exception as e:
            print("ERROR")
            print("  ->", e)

    # Run the simulation, collecting voltages
    print('Simulation', end=' ')
    model = arbor.single_cell_model(cell)
    model.probe('voltage', 'center', frequency=10000)
    model.run(tfinal=1500)
    print('DONE')
    for t in model.traces:
        ts = t.time[:]
        volts[i] = t.value[:]
        break
    spikes = np.array(model.spikes)
    count = len(spikes[(spikes >= 200) & (spikes <= 300)])
    print('Counted spikes', count)
    counts[i] = count

# Plot voltage traces
fg, ax = plt.subplots()
for k, v in volts.items():
    ax.plot(ts, v, ls='-', label=f"I={k}mA Spikes={counts[k]}")
    ax.set_xlim(left=200, right=300)
ax.set_ylabel('U/mV')
ax.set_xlabel('t/ms')
ax.legend()
plt.show()
