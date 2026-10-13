
"""
Groenheide et al. 2026
Barrovian metamorphism signals syn-orogenic delamination

Comparing the heat flux and temperature against time for varying heating windows


"""
import numpy as np
import matplotlib.pyplot as plt

#%%
###########################################################
# Model geometry & material parameters 
###########################################################

L1 = 25e3   # crystalline basement
L2 = 5e3    # sedimentary cover
L3 = 25e3   # overburden

def material_properties_standalone(z):
    Nz = len(z)
    hcond = np.zeros(Nz)
    for i in range(Nz):
        if z[i] < L1:                   # basement
            hcond[i] = 3.27
        elif z[i] < L1 + L2:            # cover
            hcond[i] = 2.12
        else:                           # overburden
            hcond[i] = 2.12
    return hcond

#%%
###########################################################
#Loading the files
###########################################################

def load_results_from_ascii(path):
    with open(path, "r") as f:
        header = f.readline().strip().split("\t")

    z_km = np.array([float(x) for x in header[1:]])
    z_m = z_km * 1000.0

    data = np.loadtxt(path, skiprows=1)
    time_s    = data[:, 0]
    results_K = data[:, 1:]

    Nz = z_m.size
    h = (z_m[-1] - z_m[0]) / (Nz - 1)
    hcond = material_properties_standalone(z_m)

    return time_s, z_m, h, hcond, results_K

#%%
###########################################################
# FLux calculation
###########################################################

def calculate_flux(time_s, results_K, z_m, h, hcond, depth_km):
    domain_depth_m  = z_m[-1]
    domain_depth_km = domain_depth_m / 1000.0

    target_depth_m   = depth_km * 1000.0
    target_height_m  = domain_depth_m - target_depth_m
    i = np.argmin(np.abs(z_m - target_height_m))
    actual_depth_km = domain_depth_km - z_m[i] / 1000.0

    if np.isclose(actual_depth_km, domain_depth_km):
        flux = hcond[0] * (results_K[:, 0] - results_K[:, 1]) / h
    else:
        if i == 0 or i == len(z_m) - 1:
            raise ValueError("Selected depth lies at a model boundary.")
        T_deeper    = results_K[:, i - 1]
        T_shallower = results_K[:, i + 1]
        k_local = hcond[i]
        flux = k_local * (T_deeper - T_shallower) / (2.0 * h)

    flux_mWm2 = flux * 1000.0
    return flux_mWm2, actual_depth_km

#%%
###########################################################
# Making tghe plot
###########################################################

def plot_flux_T_multi(ascii_files, lags_Myr, depth_km, out_png):
   

    Myr_s = 1e6 * 365.25 * 24 * 3600
    colors = ["tab:green", "tab:blue", "tab:orange", "tab:red", "tab:purple"]

    fig, ax1 = plt.subplots(dpi=300)
    ax2 = ax1.twinx()

    for path, lag, color in zip(ascii_files, lags_Myr, colors):

        time_s, z_m, h, hcond, results_K = load_results_from_ascii(path)
        t_Myr = time_s / Myr_s

        flux_mWm2, actual_depth_km = calculate_flux(
            time_s, results_K, z_m, h, hcond, depth_km
        )

        domain_depth_m   = z_m[-1]
        target_height_m  = domain_depth_m - depth_km * 1000.0
        depth_index      = np.argmin(np.abs(z_m - target_height_m))
        T_C              = results_K[:, depth_index] - 273.15

        # Temperature: solid line, labelled by lag
        ax2.plot(t_Myr, T_C, color=color, linestyle="-",
                 linewidth=2.0, label=f"{lag} Myr")

        # Flux: dashed line, no legend entry
        ax1.plot(t_Myr, flux_mWm2, color=color, linestyle="--",
                 linewidth=1.5, label="_nolegend_")
    plt.axhline(y=700, color='black', linestyle='--', linewidth=0.7)
    plt.axhline(y=550, color='red', linestyle='--', linewidth=0.7)
    ax1.set_xlabel("Time (Myr)")
    ax1.set_ylabel("Heat flux (mW/m²)")
    ax2.set_ylabel("Temperature (°C)")
    ax1.grid(True)

    # Legend: only temperature curves (lags)
    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend(h2, l2, title="Heating Window", loc="best")

    #plt.title(f"Heat Flux and Temperature Evolution at {depth_km:.0f} km\n"
             # "for Different Heating Windows")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.show()

#%%
###########################################################
# 4. Main: 
###########################################################

if __name__ == "__main__":

    ascii_files = [
        "M_700.0_1300_tl0.ascii",
        "M_700.0_1300_tl3.ascii",
        "M_700.0_1300_tl5.ascii",
        "M_700.0_1300_tl7.ascii",
        "M_700.0_1300_tl10.ascii"
    ]

    lags_Myr = [0, 3, 5, 7, 10]

    depth_km = 30.0
    out_png  = f"flux_T_multi_{depth_km:.0f}km.png"

    plot_flux_T_multi(ascii_files, lags_Myr, depth_km, out_png)