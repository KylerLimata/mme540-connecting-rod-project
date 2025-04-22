import simulate as sim

## Define piston/engine parameters
params = {
    'B': 0.08, # m
    'S': 0.096, # m
    'r': 0.2, # m
    'CR': 10,
    'T1': 21, # Celcius
    'P1': 101.325*(10**3), # Pa
    'T3': 1871 # Celcius
}

## Stress functions
funcs = []

## 
npoints = 200

results = sim.simulate_connecting_rod(params, funcs, npoints)

sim.plot_results(results)