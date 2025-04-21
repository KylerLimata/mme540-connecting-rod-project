import numpy as np

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
    theta_rod = np.zeros_like(theta_crank)

    for i in range(theta_crank.size):
        theta_cranki = theta_crank[i]
        si = s[i]
        d2 = 0
        sign = 1

        # Wrap the crank angle
        while theta_cranki < 0:
            theta_cranki += np.pi
        while theta_cranki > np.pi*2:
            theta_cranki -= np.pi
        
        if theta_cranki > np.pi:
            theta_cranki = 2*np.pi - theta_cranki
            sign = -1
        
        if theta_cranki >= 0 and theta_cranki <= np.pi/2.0:
            d2 = si - a*np.sin(np.pi/2 - theta_cranki)
        elif theta_cranki > np.pi/2.0 and theta_cranki <= np.pi:
            d2 = si + a*np.sin(theta_cranki - np.pi/2)

        theta_rod[i] = sign*(np.arcsin(d2/r) - np.pi/2)

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
    npoints : Number of points to simulate at

    Returns
    -------
    Dictionary : Computed principle stresses
    """

    output = {} # Outpput data

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
    theta_crank_compression = np.linspace(np.pi, 2*np.pi, npoints/2)
    kinematics_data_compression = piston_kinematics(B, S, r, theta_crank_compression)

    theta_crank_power = np.linspace(2*np.pi, 3*np.pi, npoints/2)
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
    output['V1'] = V1
    output['V2'] = V2
    output['V3'] = V3
    output['V4'] = V4
    output['P2'] = P2
    output['P3'] = P3
    output['P4'] = P4
    output['T2'] = T2
    output['T4'] = T4

    ## Find pressure for the compression stroke
    V_compression = kinematics_data_compression['Vd']
    P_compression = P1*(V1**k)/(V_compression**k) + Vc

    ## Find pressure for the power stroke
    V_power = kinematics_data_power['Vd'] + Vc
    P_power = (P3*V3**k)/(V_power**k)

    theta_rod_compression = kinematics_data_compression['theta_rod']
    theta_rod_power = kinematics_data_power['theta_rod']

    theta_crank = np.concat((theta_crank_compression, theta_crank_power))
    theta_rod = np.concat((theta_rod_compression, theta_rod_power))
    V = np.concat((V_compression, V_power))
    P = np.concat((P_compression, P_power))
    output['theta_crank'] = theta_crank
    output['theta_rod'] = theta_rod
    output['V'] = V
    output['P'] = P

    ## Computes load F in x and y
    A = np.pi/4*B**2 # Piston Head Area
    F = P*A # Force
    Fx = F*np.sin(theta_rod)
    Fy = F*np.cos(theta_rod)
    output['F'] = F
    output['Fx'] = Fx
    output['Fy'] = Fy

    return output