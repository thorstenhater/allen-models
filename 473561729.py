#!/usr/bin/env python3

import matplotlib.pyplot as plt

from collections import defaultdict
from pprint import pprint
import arbor
import json

# Load tree and build morphology
tree = arbor.load_swc('473561729/Gad2-IRES-Cre_Ai14_IVSCC_-172679.03.01.01_471076778_m.swc')
morph = arbor.morphology(tree, spherical_root=True)

# Label regions
labels = arbor.label_dict({'soma': '(tag 1)',
                           'axon': '(tag 2)',
                           'dend': '(tag 3)',
                           'apic': '(tag 4)',
                           'center': '(location 0 0.5)'})

# Build cell and attach Clamp and Detector
cell = arbor.cable_cell(tree, labels)

cell.place('center', arbor.spike_detector(-10))      # from tutorial
cell.place('center', arbor.iclamp(200, 1000, 0.270)) # from modeldb

# read json file and proceed to set parameters and mechanisms
with open('473561729/473561729_fit.json') as fd:
    fit = json.load(fd)


# set global values
properties = fit['conditions'][0]
T  = properties['celsius'] + 273.15
Vm = properties['v_init']
cell.set_properties(tempK=T, Vm=Vm)

# Set reversal potentials
for kv in properties['erev']:
    region = kv['section']
    for k, v in kv.items():
        if k == 'section':
            continue
        ion = k[1:]
        print(f'Setting reversal potential for species {ion} to {v}')
        cell.paint(region, arbor.ion(ion, rev_pot=float(v)))

# Setup mechanisms and parameters

## collect parameters in dict
mechs = defaultdict(dict)

### Passive parameters
pas = fit['passive'][0]
e_pas  = pas['e_pas']
ra_pas = pas['ra']

for v in pas['cm']:
    cm_pas = v['cm']
    region = v['section']
    print(f"{region:10} -> {'pas':10}: rL={ra_pas}, Vm={e_pas} cm={cm_pas}")
    cell.paint(region, cm=cm_pas, rL=ra_pas, Vm=e_pas,)

### Remaining parameters
for block in fit['genome']:
    mech   = block['mechanism'] or 'pas'
    region = block['section']
    name   = block['name'][:-(len(mech) + 1)]
    mechs[(mech, region)][name] = block['value']

## Now paint the cell using the dict
for (mech, region), vs in mechs.items():
    print(f"{region:10} -> {mech:10}: {str(vs):60}", end=' ')
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
model.probe('voltage', 'center', frequency=100)
model.run(tfinal=1200)
print('DONE')

# Plot voltage traces
fg, ax = plt.subplots()
for t in model.traces:
    print(t.time, t.value)
    ax.plot(t.time, t.value, ls='-')
    break
plt.show()
