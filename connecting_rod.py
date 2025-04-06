import numpy as np
import matplotlib.pyplot as plt

def displacement_volume(B, S, r, theta):
    """
    Calculate the displacement volume of the piston 
    at a given crank angle angle.
    
    Parameters
    ----------
    B : The piston bore
    S : The piston stroke
    r : Connecting rod length
    theta : The crank angle

    Returns
    -------
    ArrayLike : The piston displacement volume
    """
    a = S/2 # Crank Offset
    s = a*np.cos(theta) + np.sqrt(r**2 - (a**2)*(np.sin(theta)**2)) # Displacement

    return (r + a - s)*(np.pi/4)*(B**2)


theta = np.linspace(np.pi, 2*np.pi, 360)
B = 80 # Bore (mm)
BS = 1.1 # Bore/Stroke ratio
S = B/BS # Stroke (mm)
r = 200 # Connecting rod length (mm)

Vd = displacement_volume(B, S, r, theta)*(10**-6)

fig, ax = plt.subplots()

ax.plot(theta, Vd)
ax.set_xlabel("Crank Angle")
ax.set_ylabel("Displacement Volume")

# plt.show()

# print(Vd[-1])

def compute_pv(params):
    ## Clearance volume
    B = params['B'] # Bore
    S = params['S'] # Stroke
    r = params['r'] # Connecting rod length
    CR = params['CR'] # Compression Ratio
    Vc = (1/9)*(np.pi*S*B**2)/4 # Clearance volume

    ## Compute P, V, and T at all four points
    k = 1.4 # Specific heat ratio
    T1 = params['T1'] # Temperature at 1
    P1 = params['P1'] # Pressure at 1
    V1 = Vc + (np.pi*S*B**2)/4 # Volume at 1
    V2 = Vc # Volume at 2
    V3 = V2 # Volume at 3
    V4 = V1 # Volume at 4
    T3 = params['T3'] # Temperature at 3
    T2 = T1*CR**(k - 1) # Temperature at 2
    P2 = P1*CR**k
    P3 = P2*(T3/T2)
    T4 = T3/(CR**(k - 1)) # Temperature at 4
    P4 = (P1*T4)/T1 # Pressure at 4

    # Compute pressure for the compression stroke
    theta_compression = np.linspace(np.pi, 2*np.pi, 360)
    V_compression = displacement_volume(B, S, r, theta_compression) + Vc
    P_compression = P1*(V1**k)/(V_compression**k)

    # Compute pressure for the power stroke
    theta_power = np.linspace(2*np.pi, 3*np.pi, 360)
    V_power = displacement_volume(B, S, r, theta_power) + Vc
    P_power = (P3*V3**k)/(V_power**k)

    V = np.concat((V_compression, V_power))
    P = np.concat((P_compression, P_power))

    return { 
        'V': V, 
        'P': P, 
        'V1': V1, 
        'P1': P1, 
        'V2': V2, 
        'P2': P2, 
        'V3': V3, 
        'P3': P3, 
        'V4': V4, 
        'P4': P4
        }

params = {
    'B': B,
    'S': S,
    'r': r,
    'CR': 10,
    'T1': 21,
    'P1': 101.325,
    'T3': 1871
}

out = compute_pv(params)

# print(out['P'])

fig, ax = plt.subplots()

ax.plot(out['V'], out['P'])
ax.plot(out['V1'], out['P1'], 'rx')
ax.plot(out['V2'], out['P2'], 'rx')
ax.plot(out['V3'], out['P3'], 'rx')
ax.plot(out['V4'], out['P4'], 'rx')
ax.set_xlabel("Volume")
ax.set_ylabel("Pressure")

plt.show()