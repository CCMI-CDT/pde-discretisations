import numpy as np
import matplotlib.pyplot as pp
import os

# Solution smooth initial condition
def u_exact_smooth(x, t, a):
    return np.sin(2*np.pi*(x - a*t))

# Solution w/ discontinuous initial condition
def u_exact_discont(x, t, a):
    return np.where((x - a*t) % 1.0 > 0.5, 1.0, 0.0)

def residual(x, t, dt, dx, a, u_exact):
    #u_n^m, u_{n+1}^m 
    u_n   = u_exact(x, t, a)
    u_np1 = u_exact(x, t + dt, a)
    
    #u_{m+1}^{n+1}, u_{m-1}^{n+1}
    u_np1_mp1 = np.roll(u_np1, -1)
    u_np1_mm1 = np.roll(u_np1, 1)
    
    # Time discretisation
    time_term = (u_np1 - u_n) / dt 
    
    # Centered difference discretisation
    space_term = a * (u_np1_mp1 - u_np1_mm1) / (2*dx)
    
    # Residual
    tau = np.abs(time_term + space_term)
    
    return tau

def consistency_dt(u_exact):
    a = 1.0
    L = 1.0
    M = 400 # Set dx small enough to have no effect on dt error
    dx = L / M
    x = np.arange(0, L, dx)
    t = 0.0
    
    dts = np.linspace(0.001, 0.1, 100) #time discretisations
    
    errors = []
    for dt in dts:
        tau = residual(x, t, dt, dx, a, u_exact)
        err = np.linalg.norm(tau, np.inf)
        errors.append(err)
    errors = np.array(errors)
    
    # line of best fit
    log_dt = np.log(dts)
    log_err = np.log(errors)
    slope, intercept = np.polyfit(log_dt, log_err, 1)
    
    results_dir = 'consistency_results'
    os.makedirs(results_dir, exist_ok=True)
    
    pp.figure(figsize=(6,4))
    pp.loglog(dts, errors, 'b-', linewidth=2, label='Residual')
    pp.loglog(dts, np.exp(intercept)*dts**slope, 'r--', label=f'Fit slope ≈ {slope:.2f}')
    pp.xlabel(r'$\Delta t$')
    pp.ylabel('Residual (sup norm)')
    pp.legend()
    filepath = os.path.join(results_dir, "dt.png")
    pp.savefig(filepath)
    pp.close()

def consistency_dx(u_exact):
    a = 1.0
    L = 1.0
    t = 0.0
    dt = 1e-7 # Set dt small enough to have no effect on dx error
    
    Ms = np.arange(50, 2000) # Spatial discretisations
    
    errors = []
    for M in Ms:
        dx = L / M
        x = np.linspace(0, L, M, endpoint=False)
        tau = residual(x, t, dt, dx, a, u_exact)
        err = np.linalg.norm(tau, np.inf)
        errors.append(err)
    errors = np.array(errors)
    
    #line of best fit 
    log_dx = np.log(1.0 / Ms)
    log_err = np.log(errors)
    slope, intercept = np.polyfit(log_dx, log_err, 1)
    
    results_dir = 'consistency_results'
    os.makedirs(results_dir, exist_ok=True)
    
    pp.figure(figsize=(6,4))
    pp.loglog(L/Ms, errors, 'b-', linewidth=2, markersize=6, label='Residual')
    pp.loglog(L/Ms, np.exp(intercept)*(L/Ms)**slope, 'r--', label=f'Fit slope ≈ {slope:.2f}')
    pp.xlabel(r'$\Delta x$')
    pp.ylabel('Residual (sup norm)')
    pp.legend()
    filepath = os.path.join(results_dir, "dx.png")
    pp.savefig(filepath)
    pp.close()

consistency_dt(u_exact_smooth)