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

def simulate_connecting_rod(B, S, r, funcs):
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
    B : Piston Bore
    S : Piston Stroke
    r : Connecting Rod Length
    funcs : A list of functions for the normal and shear stresses at each point

    Returns
    -------
    Dictionary : Computed principle stresses
    """


    pass