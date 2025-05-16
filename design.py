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
    'd_out': 0.06, # m
    'w_beam': 0.05, # m, outer width of beam section
    't_beam': 0.03, # m, outer length of the beam
    'w_web': 0.03, # m, width of the web 
    't_web': 0.015, # m
    'w_base': 0.1, # m, width of piston base
    'r_base_fillet': 0.01, # m, fillet where the beam meets base
    'r_web_fillet': 0.0045, # m, fillet inside the beam
     # Stress Concentrations
    'kt': {
        'axial': [
            1.8, # Point 1
            1.65, # Point 2
            2.15, # Point 3
            1.0, # Point 4, not used
            1.0, # Point 5
        ],
        'bending': [
            1.5, # Point 1
            1.0, # Point 2, not used
            1.2, # Point 3
        ]
    },
    # Otto Cycle Parameters
    'T1': 21, # Celcius
    'P1': 101.325*(10**3), # Pa
    'T4': 1871 # Celcius
}

## 
npoints = 100
results = sim.simulate_new(params, npoints)

sim.plot_new(results)