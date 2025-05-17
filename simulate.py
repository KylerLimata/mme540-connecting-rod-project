import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def piston_kinematics(B, S, r, theta_crank):
    """
    Calculates the displacement volume and connecting
    rod angle for the given piston geometry at the
    passed values of theta.

    Parameters
    ----------
    B : Piston Bore
    S : Piston Stroke
    r : Connecting Rod Length
    theta_crank : The crank shaft angles

    Returns
    -------
    
    """
    ## Find the displacement volume
    a = S/2 # Crank Offset
    s = a*np.cos(theta_crank) + np.sqrt(r**2 - (a**2)*(np.sin(theta_crank)**2)) # Displacement
    Vd = (r + a - s)*(np.pi/4)*(B**2) # Displacement Volume

    ## Find the angle of the connecting rod from the y-axis
    inside = (s**2 + r**2 - a**2)/(2*s*r)
    inside_clipped = np.clip(inside, -1.0, 1.0)
    theta_rod = np.arccos(inside_clipped)

    return {'Vd': Vd, 'theta_rod': theta_rod }

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
    M = Fx*(r_rod - 0.5*w_base - r_base_fillet - 0.5*w_web)
    I = (w_web*t_web**3)/12 + ((t_beam**3)/12)*(w_beam - w_web)
    sigma_x = (M*w_beam)/(2*I)
    sigma_y = Fy/A_cross
    tau_xy = 1.5*(Fx/w_beam*t_web)

    return sigma_x, sigma_y, tau_xy

def stress_at_4(params, Fx, Fy):
    kt_axial = params['kt']['axial'][4]
    d_pin = params['d_pin']
    d_out = params['d_out']
    t_beam = params['t_beam']

    A_cross = ((d_out-d_pin)*t_beam)

    sigma_x = np.zeros_like(Fx)
    sigma_y = kt_axial*Fy/A_cross
    tau_xy = 2*Fx/A_cross

    return sigma_x, sigma_y, tau_xy
    """
    Performs a numerical simulation of the loading on a
    piston connecting rod over the compression and
    power strokes of the ideal otto cycle.

    Uses the passed functions to find the normal and
    shear stresses at any given points and computes
    the maximum principle stress at each and the
    average of the maximum principle stresses.

    2D stress conditions are assumed.

    Parameters
    ----------
    params: Piston Parameters
    funcs : A list of functions for the normal and shear stresses at each point
    npoints : Number of points to simulate at for each stroke

    Returns
    -------
    Dictionary : Computed principle stresses
    """

    results = {} # Outpput data

    ## Unpack parameters
    B = params['B'] # Bore
    S = params['S'] # Stroke
    r = params['r'] # Connecting rod length
    CR = params['CR'] # Compression Ratio
    T1 = params['T1'] # Temperature at 1
    P1 = params['P1'] # Pressure at 1
    T3 = params['T3'] # Temperature at 3
    Vc = (1/(CR - 1))*(np.pi*S*B**2)/4 # Clearance volume

    ## Calculate kinematics for compression and power stroke
    theta_crank_compression = np.linspace(np.pi, 2*np.pi, npoints)
    kinematics_data_compression = piston_kinematics(B, S, r, theta_crank_compression)

    theta_crank_power = np.linspace(2*np.pi, 3*np.pi, npoints)
    kinematics_data_power = piston_kinematics(B, S, r, theta_crank_power)

    ## Compute P, V, and T at all four points
    k = 1.4 # Specific heat ratio
    V1 = Vc + (np.pi*S*B**2)/4 # Volume at 1
    V2 = Vc # Volume at 2
    V3 = V2 # Volume at 3
    V4 = V1 # Volume at 4
    T2 = T1*CR**(k - 1) # Temperature at 2
    P2 = P1*CR**k
    P3 = P2*(T3/T2)
    T4 = T3/(CR**(k - 1)) # Temperature at 4
    P4 = (P1*T4)/T1 # Pressure at 4
    results['V1'] = V1
    results['V2'] = V2
    results['V3'] = V3
    results['V4'] = V4
    results['P1'] = P1
    results['P2'] = P2
    results['P3'] = P3
    results['P4'] = P4
    results['T2'] = T2
    results['T4'] = T4

    ## Find pressure for the compression stroke
    V_compression = kinematics_data_compression['Vd'] + Vc
    P_compression = P1*(V1**k)/(V_compression**k)

    ## Find pressure for the power stroke
    V_power = kinematics_data_power['Vd'] + Vc
    P_power = (P3*V3**k)/(V_power**k)

    theta_rod_compression = -kinematics_data_compression['theta_rod']
    theta_rod_power = kinematics_data_power['theta_rod']

    theta_crank = np.concat((theta_crank_compression, theta_crank_power))
    theta_rod = np.concat((theta_rod_compression, theta_rod_power))
    V = np.concat((V_compression, V_power))
    P = np.concat((P_compression, P_power))
    results['theta_crank'] = theta_crank
    results['theta_rod'] = theta_rod
    results['V'] = V
    results['P'] = P

    ## Compute load F in x and y
    A = np.pi/4*B**2 # Piston Head Area
    F = P*A # Force
    Fx = F*np.sin(theta_rod)
    Fy = -F*np.cos(theta_rod)
    results['F'] = F
    results['Fx'] = Fx
    results['Fy'] = Fy

    ## Evaluate each stress equation and
    ## Compute sigma 1
    funcs = [stress_at_1, stress_at_2, stress_at_3, stress_at_4]
    x_stresses = []
    y_stresses = []
    tau_stresses = []
    principal_stresses = []

    for func in funcs:
        sigma_x, sigma_y, tau_xy = func(params, Fx, Fy)
        sigma_1 = (sigma_x + sigma_y)/2 + np.sqrt(((sigma_x - sigma_y)/2)**2 + tau_xy**2)
        
        x_stresses.append(sigma_x)
        y_stresses.append(sigma_y)
        tau_stresses.append(tau_xy)
        principal_stresses.append(sigma_1)
    

    stresses = {
        'sigma_x': x_stresses,
        'sigma_y': y_stresses,
        'tau_xy': tau_stresses,
        'principal': principal_stresses
    }
    results['stresses'] = stresses

    return results

def simulate_rod(params, npoints):
    """
    Performs a numerical simulation of the loading on a
    piston connecting rod over all four strokes of the
    ideal Otto cycle.

    Uses the passed functions to find the normal and
    shear stresses at any given points and computes
    the maximum principle stress at each and the
    average of the maximum principle stresses.

    2D stress conditions are assumed.

    It is assumed that point 1 on the PV diagram is
    the start of the intake stroke, points 2-3 are the
    start and end of the compression stroke, 3-4 are
    the start and end of the power stroke, and 6-1
    are the start end end of the exhaust stroke.

    Parameters
    ----------
    params : Simulation Parameters
    npoints : Number of points to simulate at for each stroke

    Returns
    -------
    Dictionary : Computed principle stresses
    """

    results = {} # Outpput data

    ## Unpack parameters
    B = params['B'] # Bore
    S = params['S'] # Stroke
    r = params['r'] # Connecting rod length
    CR = params['CR'] # Compression Ratio
    T1 = params['T1'] # Temperature at 1
    P1 = params['P1'] # Pressure at 1
    T4 = params['T4'] # Temperature at 4
    Vs = (np.pi*S*B**2)/4 # Swept volume
    Vc = (1/(CR - 1))*Vs # Clearance volume

    ## Calculate kinematics for all four strokes
    theta_crank_intake = np.linspace(0.0, np.pi, npoints)
    kinematics_data_intake = piston_kinematics(B, S, r, theta_crank_intake)

    theta_crank_compression = np.linspace(np.pi, 2*np.pi, npoints)
    kinematics_data_compression = piston_kinematics(B, S, r, theta_crank_compression)

    theta_crank_power = np.linspace(2*np.pi, 3*np.pi, npoints)
    kinematics_data_power = piston_kinematics(B, S, r, theta_crank_power)

    theta_crank_exhaust = np.linspace(3*np.pi, 4*np.pi, npoints)
    kinematics_data_exhaust = piston_kinematics(B, S, r, theta_crank_exhaust)

    ## Compute pressure at all 6 points
    k = 1.4 # Specific heat ratio
    V1 = Vc
    # Intake Stroke, points 1-2
    V2 = Vc + Vs
    P2 = P1
    T2 = T1
    # Compression Stroke, points 2-3
    V3 = Vc
    P3 = P1*CR**k
    T3 = T1*CR**(k - 1)
    # Combustion, points 3-4
    V4 = Vc
    P4 = P3*(T4/T3)
    # Power stroke, points 4-5
    V5 = Vc + Vs
    P5 = P4/(CR**k)
    T5 = P5/(CR**(k-1))
    # Head Rejection, points 5-6
    V6 = Vc + Vs
    P6 = P2
    T6 = T2

    # Record to results
    results['V1'] = V1
    results['V2'] = V2
    results['V3'] = V3
    results['V4'] = V4
    results['V5'] = V5
    results['V6'] = V6
    results['P1'] = P1
    results['P2'] = P2
    results['P3'] = P3
    results['P4'] = P4
    results['P5'] = P5
    results['P6'] = P6
    results['T1'] = T1
    results['T2'] = T2
    results['T3'] = T3
    results['T4'] = T4
    results['T5'] = T5
    results['T6'] = T6

    ## Pressure and volume for intake stroke (1-2)
    V_intake = np.linspace(V1, V2, len(theta_crank_intake))
    P_intake = np.ones_like(theta_crank_intake)*P1

    ## Pressure and volume for compression stroke (2-3)
    V_compression = kinematics_data_compression['Vd'] + Vc
    P_compression = P2*(V2**k)/(V_compression**k)

    ## Pressure and volume for power stroke (4-5)
    V_power = kinematics_data_power['Vd'] + Vc
    P_power = (P4*V4**k)/(V_power**k)

    ## Pressure and valume for exhaust stroke (6-1)
    V_exhaust = np.linspace(V6, V1, len(theta_crank_exhaust))
    P_exhaust = np.ones_like(theta_crank_intake)*P6

    ## Rod angle
    theta_rod_intake = kinematics_data_intake['theta_rod']
    theta_rod_compression = -kinematics_data_compression['theta_rod']
    theta_rod_power = kinematics_data_power['theta_rod']
    theta_rod_exhaust = -kinematics_data_exhaust['theta_rod']

    theta_crank = np.concat((
        theta_crank_intake, 
        theta_crank_compression, 
        theta_crank_power,
        theta_crank_exhaust
        ))
    theta_rod = np.concat((
        theta_rod_intake, 
        theta_rod_compression, 
        theta_rod_power,
        theta_rod_exhaust
        ))
    V = np.concat((
        V_intake,
        V_compression,
        V_power,
        V_exhaust
    ))
    P = np.concat((
        P_intake,
        P_compression,
        P_power,
        P_exhaust
    ))

    results['theta_crank'] = theta_crank
    results['theta_rod'] = theta_rod
    results['P'] = P
    results['V'] = V

    ## Compute load F in x and y
    A = np.pi/4*B**2 # Piston Head Area
    F = P*A # Force magnitude
    Fx = F*np.sin(theta_rod)
    Fy = -F*np.cos(theta_rod)
    results['F'] = F
    results['Fx'] = Fx
    results['Fy'] = Fy

    ## Evaluate each stress equation and
    ## Compute sigma 1
    funcs = [stress_at_1, stress_at_2, stress_at_3, stress_at_4]
    x_stresses = []
    y_stresses = []
    tau_stresses = []
    principal_stresses = []

    for func in funcs:
        sigma_x, sigma_y, tau_xy = func(params, Fx, Fy)
        sigma_1 = (sigma_x + sigma_y)/2 + np.sqrt(((sigma_x - sigma_y)/2)**2 + tau_xy**2)
        
        x_stresses.append(sigma_x)
        y_stresses.append(sigma_y)
        tau_stresses.append(tau_xy)
        principal_stresses.append(sigma_1)
    

    stresses = {
        'sigma_x': x_stresses,
        'sigma_y': y_stresses,
        'tau_xy': tau_stresses,
        'sigma_1': principal_stresses
    }
    results['stresses'] = stresses

    return results
    theta_crank = results['theta_crank']

    ## Plot the PV diagram
    P = np.append(results['P'], results['P1'])
    V = np.append(results['V'], results['V1'])
    fig, ax = plt.subplots()

    ax.plot(V*10**9, P*10**-3)
    ax.set_xlabel("$Volume (mm^3)$")
    ax.set_ylabel("Pressure (kPa)")
    ax.set_title("PV Diagram")

    ## Plot rod angle
    theta_crank = results['theta_crank']
    theta_rod = results['theta_rod']
    fig, ax = plt.subplots()

    ax.plot(theta_crank, theta_rod)
    ax.set_xlabel(r"$\theta_{crank} (rad)$")
    ax.set_ylabel(r"$\theta_{rod}$ (rad)")
    ax.set_title("Rod Angle vs. Crank Angle")

    ## Plot Forces
    F = results['F']*10**-3
    Fx = results['Fx']*10**-3
    Fy = results['Fy']*10**-3
    fig, ax = plt.subplots()

    ax.plot(theta_crank, F, label = 'F')
    ax.plot(theta_crank, Fx, label = 'Fx')
    ax.plot(theta_crank, Fy, label = 'Fy')
    ax.set_xlabel(r"$\theta_{crank} (rad)$")
    ax.set_ylabel("Force (kN)")
    ax.legend()

    ## Plot x Stresses
    x_stresses = results['stresses']['sigma_x']
    fig, ax = plt.subplots()
    i = 0

    for sigma_x in x_stresses:
        i = i + 1

        ax.plot(theta_crank, sigma_x, label=r"$\sigma_{n,x}$".replace('n', f'{i}'))

    ax.set_xlabel(r"$\theta_{crank} (rad)$")
    ax.set_ylabel("Stress")
    ax.set_title(r"$\sigma_x$")
    ax.legend()

    ## Plot y Stresses
    y_stresses = results['stresses']['sigma_y']
    fig, ax = plt.subplots()
    i = 0

    for sigma_y in y_stresses:
        i = i + 1

        ax.plot(theta_crank, sigma_y, label=r"$\sigma_{n,y}$".replace('n', f'{i}'))

    ax.set_xlabel(r"$\theta_{crank} (rad)$")
    ax.set_ylabel("Stress")
    ax.set_title(r"$\sigma_y$")
    ax.legend()

    ## Plot Sheat Stresses
    tau_stresses = results['stresses']['tau_xy']
    fig, ax = plt.subplots()
    i = 0

    for tau_xy in tau_stresses:
        i = i + 1

        ax.plot(theta_crank, tau_xy, label=r"$\tau_{n,xy}$".replace('n', f'{i}'))

    ax.set_xlabel(r"$\theta_{crank} (rad)$")
    ax.set_ylabel("Stress")
    ax.set_title(r"$\tau_xy$")
    ax.legend()

    ## Plot Principal Stresses
    principal_stresses = results['stresses']['principal']
    fig, ax = plt.subplots()
    i = 0

    for sigma_x in principal_stresses:
        i = i + 1

        ax.plot(theta_crank, sigma_x, label=r"$\sigma_{n,1}$".replace('n', f'{i}'))

    ax.set_xlabel(r"$\theta_{crank} (rad)$")
    ax.set_ylabel("Stress")
    ax.set_title("Principal Stresses")
    ax.legend()

    plt.show()

def plot_results(results):
    """
    Plots the simulation results for immediate viewing.
    Produces a plot for the PV diagram, a plot for the
    load vs. crankshaft angle, and a plot of all the
    stresses.

    Parameters
    ----------
    results : connecting rod simulation results
    """
    theta_crank = results['theta_crank']

    ## Plot the PV diagram
    P = results['P']
    V = results['V']
    fig, ax = plt.subplots()

    ax.plot(V*10**9, P*10**-6)
    ax.set_xlabel("$Volume (mm^3)$")
    ax.set_ylabel("Pressure (MPa)")
    ax.set_title("PV Diagram")

    ## Plot the loads
    F = results['F']
    Fx = results['Fx']
    Fy = results['Fy']
    fig, ax = plt.subplots()

    ax.plot(theta_crank, F*10**-3, label = 'F')
    ax.plot(theta_crank, Fx*10**-3, label = 'Fx')
    ax.plot(theta_crank, Fy*10**-3, label = 'Fy')
    ax.set_xlabel(r"$\theta_{crank} (rad)$")
    ax.set_ylabel("Force (kN)")
    ax.legend()

    ## Plot the stresses
    titles = [r"$\sigma_x$", r"$\sigma_y$", r"$\tau_{xy}$", r"$\sigma_1$"]
    labels = ['Point 1', 'Point 2', 'Point 3', 'Point 4']
    colors = plt.cm.viridis(np.linspace(0, 1, 4))
    sigma_x = results['stresses']['sigma_x']
    sigma_y = results['stresses']['sigma_y']
    tau_xy = results['stresses']['tau_xy']
    sigma_1 = results['stresses']['sigma_1']
    data = [sigma_x, sigma_y, tau_xy, sigma_1]

    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)

    for i, ax in enumerate(axs.flat):
        for j, color in enumerate(colors):
            y = data[i][j]
            ax.plot(theta_crank, y, label=labels[j], color=color)
        ax.set_title(titles[i])
    
    
    # Extract legend handles and labels from the last axis (all should be the same)
    handles, labels = axs[0, 0].get_legend_handles_labels()

    # Add a single legend to the figure
    fig.legend(handles, labels, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.01))

    # Adjust layout so the legend doesn't overlap
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)

    plt.show()

def save_results(results, name):
    """
    Saves the simulation results to several
    different files. Each filename and the
    directory the files are saved in is
    suffixed with the passed name parameter

    "loads_{name}.csv": The x and y loads over
    timesteps for easy copy-pasting into
    an ANSYS transient structural simulation

    "stresses_{name}.csv": The computed stress data

    "report_{name}.txt": A report file detailing the
    max and average principle stress at each point

    Parameters
    ----------
    results : connecting rod simulation results
    name : the name of the design
    """
    ## Create the results folder if it doesn't already exist
    if not os.path.exists(f"results_{name}"):
        os.mkdir(f"results_{name}")

    ## Save the loads for copy-pasting into ANSYS Transient
    Fx = results['Fx']
    Fy = results['Fy']
    data = {
        'step': np.linspace(1, len(Fx), len(Fx)),
        'Fx': Fx,
        'Fy': Fy
    }
    df = pd.DataFrame(data)
    df.to_csv(f"results_{name}/loads_{name}.csv", index=False)

    ## Save the stresses
    sigma_x = results['stresses']['sigma_x']
    sigma_y = results['stresses']['sigma_y']
    tau_xy = results['stresses']['tau_xy']
    sigma_1 = results['stresses']['sigma_1']
    data = {
        'thetacrank': results['theta_crank'],
        'Fx': Fx,
        'Fy': Fy
    }

    for i in range(len(sigma_1)):
        data[f"sigmax{i + 1}"] = sigma_x[i]
        data[f"sigmay{i + 1}"] = sigma_y[i]
        data[f"tauxy{i + 1}"] = tau_xy[i]
        data[f"sigma1{i + 1}"] = sigma_1[i]

    df = pd.DataFrame(data)

    df.to_csv(f"results_{name}/stresses_{name}.csv", index=False)

    ## Write report
    with open(f"results_{name}/report_{name}.txt", "w") as file:
        for i, sigma_1i in enumerate(sigma_1):
            avg_sigma_1i = np.mean(sigma_1i)*10**-6
            max_sigma_1i = np.max(sigma_1i)*10**-6
            file.write(f"Point {i+1}:\n  mean stress = {avg_sigma_1i:.3f} MPa,\n  max stress = {max_sigma_1i:.3f} Mpa\n")