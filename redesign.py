import simulate as sim
import numpy as np

## Define piston/engine parameters
params = {
    # Geometry Parameters
    'B': 0.08, # m, bore
    'S': 0.096, # m, stroke
    'r': 0.2, # m, connecting rod length
    'CR': 10, # compression ratio
    'd_pin': 0.045, # m, pin diameter
    'd_ring': 0.06, # m
    't_ring': 0.06,  # m
    'w_beam': 0.05, # m, outer width of beam section
    't_beam': 0.03, # m, outer length of the beam
    'w_web': 0.03, # m, width of the web 
    'w_base': 0.1, # m, width of piston base
    'r_base_fillet': 0.01, # m, fillet where the beam meets base
     # Stress Concentrations
    'kt': {
        'axial': [
            1.8, # Point 1
            1.3, # Point 2
            1.0, # Point 3
            1.0, # Point 4, not used
        ],
        'bending': [
            1.5, # Point 1
            1.3, # Point 2, not used
            1.0, # Point 3
            1.0, # Point 4
        ]
    },
    # Otto Cycle Parameters
    'T1': 21, # Celcius
    'P1': 101.325*(10**3), # Pa
    'T4': 1327 # Celcius
}

## Compute dimensions with respect to chosen kt
tbeam_over_tweb = 1.05
rwebfillet_over_tweb = 0.3
params['t_web'] = params['t_beam']/tbeam_over_tweb
params['r_web_fillet'] = params['t_web']*rwebfillet_over_tweb

## 
npoints = 100
results = sim.simulate_rod(params, npoints)

sim.save_results(results, "redesign")
sim.plot_results(results)