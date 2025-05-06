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
    M = Fx*(r_rod - 0.5*w_base)
    sigma_x = kt_bending*((6*M)/(t_beam*w_beam**2))
    sigma_y = kt_axial*(Fy/A_cross)
    tau_xy = 1.5*(Fx/A_cross)

    return sigma_x, sigma_y, tau_xy

def stress_at_2(params, Fx, Fy):
    kt_axial = params['kt']['axial'][1]
    t_web = params['t_web']
    w_beam = params['w_beam']

    A_cross = t_web*w_beam
    sigma_x = np.zeros_like(Fx)
    sigma_y = kt_axial*Fy/A_cross
    tau_xy = 1.5*Fx/A_cross

    return sigma_x, sigma_y, tau_xy

def stress_at_2_revised(params, Fx, Fy):
    kt_axial = params['kt']['axial'][1]
    t_web = params['t_web']
    w_web = params['w_web']
    t_beam = params['w_beam']
    w_beam = params['w_beam']

    A_cross = t_web*w_beam
    A_web = w_web*t_web
    A_total = w_beam*t_beam - w_web*(t_beam - t_web)
    A_remaining = A_total - A_web
    concentrated_percent = A_web/A_total

    sigma_x = np.zeros_like(Fx)
    sigma_y = (kt_axial*concentrated_percent*Fy)/A_web
    tau_xy = 1.5*Fx/A_cross

    return sigma_x, sigma_y, tau_xy

def stress_at_3(params, Fx, Fy):
    kt_axial = params['kt']['axial'][2]
    kt_bending = params['kt']['bending'][2]
    w_beam = params['w_beam']
    t_beam = params['t_beam']
    w_web = params['w_web']
    t_web = params['t_web']
    r_rod = params['r']
    w_base = params['w_base']
    r_base_fillet = params['r_base_fillet']

    t_hole = t_beam - t_web
    A_total = w_beam*t_beam - 2*w_web*(t_beam - t_web)
    A_hole = (w_beam - w_web)*(t_beam - t_web)
    percent_depth = A_hole/A_total
    M = (r_rod - w_base/2 - r_base_fillet - w_web/2)*Fx
    I_mid = (t_web*w_beam**3)/12
    c = w_beam/2
    
    sigma_x = kt_bending*((percent_depth*6*M)/((w_beam - w_web)*t_hole**2)) \
                + ((1 - percent_depth)*M*c)/I_mid
    sigma_y = kt_axial*(percent_depth*Fy)/((w_beam - w_web)*t_hole) \
                + ((1 - percent_depth)*Fy)/(w_beam - t_web)
    tau_xy = (1.5*Fx)/(w_beam*t_web)

    return sigma_x, sigma_y, tau_xy

def stress_at_3_revised(params, Fx, Fy):
    kt_axial = params['kt']['axial'][2]
    kt_bending = params['kt']['bending'][2]
    w_beam = params['w_beam']
    t_beam = params['t_beam']
    w_web = params['w_web']
    t_web = params['t_web']
    r_rod = params['r']
    w_base = params['w_base']
    r_base_fillet = params['r_base_fillet']
    r_web_fillet = params['r_web_fillet']

    A_fillet = (1 - 0.25*np.pi)*r_web_fillet**2
    A_total = w_beam*t_beam - w_web*(t_beam - t_web) + 4*A_fillet
    A_hole = (w_beam - w_web)*(t_beam - t_web)
    concentrated_percent = A_hole/A_total

    M = (r_rod - w_base/2 - r_base_fillet - w_web/2)*Fx
    I_mid = (t_web*w_beam**3)/12
    
    sigma_x = kt_bending*concentrated_percent*((M*w_beam)/(2*I_mid))
    sigma_y = kt_axial*(concentrated_percent*Fy)/(A_hole)
    tau_xy = (1.5*Fx)/(w_beam*t_web)

    return sigma_x, sigma_y, tau_xy


def stress_at_4(params, Fx, Fy):
    w_beam = params['w_beam']
    t_beam = params['t_beam']
    w_web = params['w_web']
    t_web = params['t_web']
    r_rod = params['r']
    w_base = params['w_base']
    r_base_fillet = params['r_base_fillet']
    r_web_fillet = params['r_web_fillet']

    A_fillet = (1 - 0.25*np.pi)*r_web_fillet**2
    
    A_cross = w_beam*t_beam - w_web*(t_beam - t_web) + 4*A_fillet
    M = Fx*(r_rod - 0.5*w_base - r_base_fillet - 0.5*w_web - 0.01)
    I = (w_web*t_web**3)/12 + ((t_beam**3)/12)*(w_beam - w_web)
    sigma_x = (M*w_beam)/(2*I)
    sigma_y = Fy/A_cross
    tau_xy = 1.5*(Fx/w_beam*t_web)

    return sigma_x, sigma_y, tau_xy

funcs = [stress_at_1, stress_at_2_revised, stress_at_3_revised, stress_at_4]

## 
npoints = 200

results = sim.simulate_connecting_rod(params, funcs, npoints)

sim.analyze_results(results)
sim.save_results(results, "final")
sim.plot_results(results)