import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

def simulate_connecting_rod(params, funcs, npoints):
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
    x_stresses = []
    y_stresses = []
    tau_stresses = []
    principal_stresses = []

    for func in funcs:
        sigma_x, sigma_y, tau_xy = func(params, Fx, Fy)
        sigma_1 = (sigma_x + sigma_y)/2 + np.sqrt(((sigma_x + sigma_y)/2)**2 + tau_xy**2)
        
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

def plot_results(results):
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

def analyze_results(results):
    x_stresses = results['stresses']['sigma_x']
    y_stresses = results['stresses']['sigma_y']
    tau_stresses = results['stresses']['tau_xy']
    principal_stresses = results['stresses']['principal']
    i = 0

    np.set_printoptions(precision=4)

    print("Maximum and Average Stresses:")
    for sigma_1 in principal_stresses:
        i = i + 1
        avg_sigma_1 = np.mean(sigma_1)*10**-6
        max_sigma_1 = np.max(sigma_1)*10**-6

        print(f"At point {i}, mean stress = {avg_sigma_1:.3f} MPa, max stress = {max_sigma_1:.3f} Mpa")

    print("Stresses at the Start of Combustion:")


def save_results(results, filename):
    row_names = ['theta', 'Fx', 'Fy']
    data = np.array([results['theta_crank'],results['Fx'],results['Fy']])
    x_stresses = results['stresses']['sigma_x']

    i = 0
    for sigma_x in x_stresses:
        i = i + 1
        data = np.concatenate((data, sigma_x), axis=1)
        row_names.append(f"point {i}")

    df = pd.DataFrame(data, index=row_names)

    df.to_csv(f"results_{filename}_sigmax.csv", index=True, header=False)
