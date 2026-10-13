# -*- coding: utf-8 -*-
"""
Groenheide et al. 2026
Barrovian metamorphism signals syn-orogenic delamination

Comparing Boundary conditions for a given tmeperature. 


"""



import numpy as np
import matplotlib.pyplot as plt
import os


#%%
###########################################################
# Experiment selection and control parameter input.
###########################################################

data_folder = os.path.dirname(os.path.abspath(__file__))  # <- same folder as this script

experiment    = 2       # 1 = IR , 2 = Rad
T_moho_choice = 700            # Take a Moho T out of model range 

#Selecting the age and depth of the plots
plot_age    = 30               # Myr  
plot_depth  = 30              # km   

# grid geometry 
L1, L2, L3 = 25e3, 5e3, 25e3   # m

models = ["CTBC", "HF_q0BC", "HF_qcBC", "HF_qLAB"]


exp_tag = "IR" if experiment == 1 else "Rad"
year = 365.25 * 24 * 3600
Myr = 1e6 * year
L_total_km = (L1 + L2 + L3) / 1000

#%%
###########################################################
# Load files 
###########################################################
def load_ascii(path):
    """
    Load one .ascii file.
    Returns:
        time  : (nt,) array, seconds
        z_km  : (nz,) array, km, ordered base(0) -> top(L)
        T     : (nt, nz) array, Kelvin
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    z_km = np.array(header[1:], dtype=float)

    data = np.loadtxt(path, skiprows=1)
    time = data[:, 0]
    T = data[:, 1:]
    return time, z_km, T


def find_file(folder, model, exp_tag, Tmoho):
    fname = f"{model}_{exp_tag}_M_{Tmoho:.1f}.ascii"
    path = os.path.join(folder, fname)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Could not find:\n  {path}\n"
            f"Expected naming pattern: {{model}}_{{exp_tag}}_M_{{Tmoho}}.ascii\n"
            f"(e.g. {model}_{exp_tag}_M_700.0.ascii)"
        )
    return path


def load_all(data_folder, models, exp_tag, T_moho_choice):
    datasets = {}
    for model in models:
        path = find_file(data_folder, model, exp_tag, T_moho_choice)
        time, z_km, T = load_ascii(path)
        datasets[model] = {"time": time, "z_km": z_km, "T": T}
    return datasets


###########################################################
# Plotting moduules
###########################################################

def plot_Tz(datasets, age_myr, T_moho_choice, exp_tag):
    plt.figure(dpi=250)
    age_s = age_myr * Myr

    for label, d in datasets.items():
        t_index = np.argmin(np.abs(d["time"] - age_s))
        T_profile = d["T"][t_index] - 273.15
        plt.plot(d["z_km"], np.flip(T_profile), label=label)

    plt.xlabel("Depth (km)")
    plt.ylabel("Temperature (\u00b0C)")
    plt.gca().invert_xaxis()
    plt.axhline(y=700, color="black", linestyle="--", linewidth=0.7)
    plt.axhline(y=550, color="red", linestyle="--", linewidth=0.7)
    plt.grid(True)
    #plt.legend()
    #plt.title(f"Temperature Profile for a Radioactive Heating Model at {age_myr} Myr")
    
    plt.savefig(f"Tz_compare_{exp_tag}_{T_moho_choice}C_{age_myr}Myr.png")
    plt.show()



def plot_Tt(datasets, depth_km, T_moho_choice, exp_tag):
    plt.figure(dpi=250)

    for label, d in datasets.items():
        depth_index = np.argmin(np.abs((d["z_km"][-1] - d["z_km"]) - depth_km))
        T_time = d["T"][:, depth_index] - 273.15
        plt.plot(d["time"] / Myr, T_time, label=label)

    plt.xlabel("Time (Myr)")
    plt.ylabel("Temperature (\u00b0C)")  # Standardized LaTeX notation
    plt.axhline(y=700, color="black", linestyle="--", linewidth=0.7)
    plt.axhline(y=550, color="red", linestyle="--", linewidth=0.7)
    plt.grid(True)
    #plt.legend()
    
    # 1. Break the long title into two lines so it fits the figure width
    #plt.title(f"Temperature Evolution of a Radioactive Heating Model\nat {depth_km} km depth")
    #
    # 2. Added bbox_inches="tight" to dynamically adjust margins and prevent cutoff
    plt.savefig(f"Tt_compare_{exp_tag}_{T_moho_choice}C_{depth_km}km.png", bbox_inches="tight")
    plt.show()


# ============================================================
# Main
# ============================================================

def main():
    datasets = load_all(data_folder, models, exp_tag, T_moho_choice)

    plot_Tz(datasets, plot_age, T_moho_choice,  exp_tag)
    plot_Tt(datasets, plot_depth, T_moho_choice,  exp_tag)


if __name__ == "__main__":
    main()