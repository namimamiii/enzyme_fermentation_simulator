import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ==========================================
# 1. ENZYME KINETICS
# ==========================================
def calculate_enzyme_speed(substrate_conc, vmax=100.0, km=5.0):
    """
    Michaelis-Menten Equation: v = (Vmax * S) / (Km + S)
    vmax: Maximum reaction velocity
    km: Substrate concentration at 1/2 Vmax
    """
    return (vmax * substrate_conc) / (km + substrate_conc)


# ==========================================
# 2. FERMENTATION KINETICS (MONOD MODEL)
# ==========================================
def monod_growth_model(time, state, mu_max, ks, yield_factor):
    """
    Calculates cell growth and nutrient depletion over time.
    state: [cells (g/L), food (g/L)]
    """
    cells, food = state

    # Stop growth when nutrients are exhausted
    if food <= 0:
        mu = 0.0
        food = 0.0
    else:
        # Monod Equation: mu = mu_max * (S / (Ks + S))
        mu = mu_max * (food / (ks + food))

    d_cells = mu * cells                  # Biomass increase rate
    d_food = -(d_cells / yield_factor)    # Nutrient consumption rate

    return [d_cells, d_food]


# ==========================================
# 3. MAIN SIMULATION RUNNER
# ==========================================
def run_simulation():
    # --- PART A: ENZYME COMPARISON (37°C vs 15°C) ---
    substrate_levels = np.linspace(0, 50, 100) # 0 to 50 mM substrate
    
    # Normal Temperature (37°C): Fast reaction speed
    v_warm = calculate_enzyme_speed(substrate_levels, vmax=150.0, km=5.0)
    
    # Cold Temperature (15°C): Slowed down enzyme activity
    v_cold = calculate_enzyme_speed(substrate_levels, vmax=40.0, km=12.0)


    # --- PART B: ORGANISM COMPARISON (E. coli vs Yeast) ---
    initial_conditions = [0.1, 20.0]  # [0.1 g/L initial cells, 20.0 g/L food]
    time_points = np.linspace(0, 15, 100) # Simulate 0 to 15 hours

    # Organism 1: E. coli (Fast Growth)
    # mu_max = 0.9 hr^-1, Ks = 0.2 g/L, Yield = 0.5
    res_ecoli = solve_ivp(
        monod_growth_model, (0, 15), initial_conditions, 
        args=(0.9, 0.2, 0.5), t_eval=time_points
    )

    # Organism 2: Baker's Yeast (Slower Growth)
    # mu_max = 0.3 hr^-1, Ks = 0.5 g/L, Yield = 0.4
    res_yeast = solve_ivp(
        monod_growth_model, (0, 15), initial_conditions, 
        args=(0.3, 0.5, 0.4), t_eval=time_points
    )


    # --- PART C: PLOTTING THE RESULTS ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Temperature Effect on Enzyme Velocity
    ax1.plot(substrate_levels, v_warm, label="Warm Temp (37°C) - High Vmax", color="red", linewidth=2)
    ax1.plot(substrate_levels, v_cold, label="Cold Temp (15°C) - Low Vmax", color="blue", linestyle="--", linewidth=2)
    ax1.set_title("Enzyme Reaction Speed vs Temperature")
    ax1.set_xlabel("Substrate Concentration (mM)")
    ax1.set_ylabel("Reaction Velocity (uM/min)")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Organism Growth Comparison
    ax2.plot(res_ecoli.t, res_ecoli.y[0], label="E. coli Biomass (Fast)", color="green", linewidth=2)
    ax2.plot(res_yeast.t, res_yeast.y[0], label="Yeast Biomass (Slow)", color="purple", linewidth=2)
    ax2.plot(res_ecoli.t, res_ecoli.y[1], label="E. coli Nutrients", color="green", linestyle=":", alpha=0.6)
    ax2.plot(res_yeast.t, res_yeast.y[1], label="Yeast Nutrients", color="purple", linestyle=":", alpha=0.6)
    
    ax2.set_title("Fermentation Growth Profile (E. coli vs Yeast)")
    ax2.set_xlabel("Time (Hours)")
    ax2.set_ylabel("Concentration (g/L)")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("biotech_comparison_plots.png", dpi=300)
    print("Success! Generated 'biotech_comparison_plots.png'")
    plt.show()

# Run script
if __name__ == "__main__":
    run_simulation()