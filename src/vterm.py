from scipy.integrate import solve_ivp
import numpy as np
from math import sqrt
from matplotlib import pyplot as plt

# Default parameters
DEFAULT_M = 1.0
DEFAULT_AIR_K = 0.12 # From projScPY2.py

def hit_ground(t, y):
    return y[2]
hit_ground.terminal = True
hit_ground.direction = -1

def func(t, y, m, air_k):
    g = 9.81
    v = sqrt(y[1]*y[1] + y[3]*y[3])
    
    # y[0] = x, y[1] = vx, y[2] = y, y[3] = vy
    f0 = y[1]                         # dx/dt = vx
    f1 = -air_k * v * y[1] / m        # dvx/dt = -k/m * v * vx
    f2 = y[3]                         # dy/dt = vy
    f3 = -air_k * v * y[3] / m - g    # dvy/dt = -k/m * v * vy - g
    
    return [f0, f1, f2, f3]

def run_simulation(m, air_k, y0, t_span):
    # Wrapper for solve_ivp to pass args
    return solve_ivp(lambda t, y: func(t, y, m, air_k), t_span, y0, 
                     events=[hit_ground], rtol=1e-8, atol=1e-8)

def check_energy_conservation():
    print("--- Checking Energy Conservation (No Air Resistance) ---")
    m = 1.0
    air_k = 0.0
    y0 = [0, 10, 0, 10]
    t_span = [0, 5]
    
    sol = run_simulation(m, air_k, y0, t_span)
    
    if len(sol.t_events[0]) > 0:
        print(f"Hit ground at t = {sol.t_events[0][0]:.3f} s")
        
    yf = sol.y
    vx = yf[1]
    y_pos = yf[2]
    vy = yf[3]
    g = 9.81
    
    ke = 0.5 * m * (vx*vx + vy*vy)
    pe = m * g * y_pos
    total_energy = ke + pe
    
    e0 = total_energy[0]
    ef = total_energy[-1]
    print(f"Initial Energy: {e0:.4f} J")
    print(f"Final Energy: {ef:.4f} J")
    print(f"Max Deviation: {np.max(np.abs(total_energy - e0)):.4e} J")
    
    plt.figure()
    plt.plot(sol.t, total_energy, label='Total Energy')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [J]')
    plt.title('Energy Conservation (Air K = 0)')
    plt.legend()
    plt.savefig('energy_conservation.png')
    print("Saved energy_conservation.png")
    print("----------------------------------------------------\n")

def find_terminal_velocity(m, air_k):
    # Drop from high altitude to reach terminal velocity
    # v_t occurs when mg = k * v^2 -> v = sqrt(mg/k)
    # We simulate to verify.
    y0 = [0, 0, 1000, 0] # Drop from 1000m
    t_span = [0, 100]
    
    sol = run_simulation(m, air_k, y0, t_span)
    
    # Terminal velocity should be the final velocity (vy)
    vy_final = sol.y[3][-1]
    return abs(vy_final)

def study_vt_vs_mass():
    print("--- Studying Terminal Velocity vs Mass ---")
    air_k = DEFAULT_AIR_K
    g = 9.81
    
    masses = np.logspace(-3, 1, 20) # 1g (1e-3) to 10kg (1e1)
    vts_sim = []
    vts_analytic = []
    
    print(f"{'Mass [kg]':<12} {'Vt Sim [m/s]':<15} {'Vt Calc [m/s]':<15} {'Diff [%]':<10}")
    
    for m in masses:
        vt_sim = find_terminal_velocity(m, air_k)
        vt_calc = sqrt(m * g / air_k)
        
        vts_sim.append(vt_sim)
        vts_analytic.append(vt_calc)
        
        diff_pct = abs(vt_sim - vt_calc) / vt_calc * 100
        print(f"{m:<12.4f} {vt_sim:<15.4f} {vt_calc:<15.4f} {diff_pct:<10.2f}")
        
    plt.figure()
    plt.loglog(masses, vts_sim, 'o', label='Simulated')
    plt.loglog(masses, vts_analytic, '-', label='Analytical (sqrt(mg/k))')
    plt.xlabel('Mass [kg]')
    plt.ylabel('Terminal Velocity [m/s]')
    plt.title('Terminal Velocity vs Mass')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    plt.savefig('vterm_plot.pdf')
    print("Saved vterm_plot.pdf")
    print("------------------------------------------\n")

if __name__ == "__main__":
    check_energy_conservation()
    study_vt_vs_mass()
