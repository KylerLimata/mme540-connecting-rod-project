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
    't_ring': 0.06,  # m
    't_beam': 0.03, # m, outer length of the beam
    'w_web': 0.03, # m, width of the web 
    'w_base': 0.1, # m, width of piston base
     # Stress Concentrations
    'kt': {
        'axial': [
            1.65, # Point 1
            1.3, # Point 2
            1.0, # Point 3, not used
            2.4, # Point 4
        ],
        'bending': [
            1.35, # Point 1
            1.0, # Point 2, not used
            1.0, # Point 3, not used
            1.2, # Point 4
        ]
    },
    # Otto Cycle Parameters
    'T1': 21, # Celcius
    'P1': 101.325*(10**3), # Pa
    'T4': 1327, # Celcius
    # SF params
    'Cm': 0.8,
    'Cst': 0.8,
    'Cr': 0.75,
    'Sn': 1160*10**6,
    'Su': 400*10**6
}

## Compute certain dimensions to minimize stress concentration
tbeam_over_tweb = 1.05
rwebfillet_over_tweb = 0.3
params['t_web'] = params['t_beam']/tbeam_over_tweb # m, thickness of the web
params['r_web_fillet'] = params['t_web']*rwebfillet_over_tweb # m, fillet inside the web

wbase_over_wbeam = 2
rbasefillet_over_wbeam = 0.3
params['w_beam'] = params['w_base']/wbase_over_wbeam # m, width of the beam section
params['r_base_fillet'] = params['w_beam']*rbasefillet_over_wbeam # m, fillet where the beam meets base

dpin_over_dring = 0.6
params['d_ring'] = params['d_pin']/dpin_over_dring # Outer diameter of the ring around the pin

print(f"t_web = {params['t_web']}")
print(f"r_web_fillet = {params['r_web_fillet']}")
print(f"w_beam = {params['w_beam']}")
print(f"r_base_fillet = {params['r_base_fillet']}")
print(f"d_ring = {params['d_ring']}")

## 
npoints = 100
results = sim.simulate_rod(params, npoints)

sim.compute_safety_factors(params, results)
sim.save_results(results, "redesign")
sim.plot_results(results)