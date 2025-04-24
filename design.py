import simulate as sim

## Define piston/engine parameters
params = {
    # Geometry Parameters
    'B': 0.08, # m
    'S': 0.096, # m
    'r': 0.2, # m
    'CR': 10,
    'w_beam': 0.05, # m
    't_beam': 0.03, # m
    't_center': 0.02, # m
    'w_base': 0.1, # m
    'r_base_fillet': 0.01, # m
     # Stress Concentrations
    'kt': {
        'axial': [
            1.175, # Points 1 and 2
        ],
        'bending': [
            1.15, # Points 1 and 2
        ]
    },
    # Otto Cycle Parameters
    'T1': 21, # Celcius
    'P1': 101.325*(10**3), # Pa
    'T3': 1871 # Celcius
}

## Stress functions
def stress_at_1(params, Fx, Fy):
    kt_axial = params['kt']['axial'][0]
    kt_bending = params['kt']['bending'][0]
    w_beam = params['w_beam']
    t_beam = params['t_beam']
    r_rod = params['r']
    w_base = params['w_base']

    A_cross = w_beam*t_beam
    M = Fx*(r_rod - w_base)
    sigma_x = 1.5*(Fx/A_cross)
    sigma_y = kt_axial*(Fy/A_cross)
    tau_xy = kt_bending*((6*M)/(t_beam*w_beam**2))

    return sigma_x, sigma_y, tau_xy

def stress_at_2(params, Fx, Fy):
    kt_axial = params['kt']['axial'][0]
    kt_bending = params['kt']['bending'][0]
    w_beam = params['w_beam']
    t_beam = params['t_beam']
    r_rod = params['r']
    w_base = params['w_base']

    A_cross = w_beam*t_beam
    M = Fx*(r_rod - w_base)
    sigma_x = 1.5*(Fx/A_cross)
    sigma_y = kt_axial*(Fy/A_cross)
    tau_xy = kt_bending*((-6*M)/(t_beam*w_beam**2))

    return sigma_x, sigma_y, tau_xy

funcs = [stress_at_1, stress_at_2]

## 
npoints = 200

results = sim.simulate_connecting_rod(params, funcs, npoints)

sim.analyze_results(results)
sim.plot_results(results)