# -*- coding: utf-8 -*-

"""""""""""""""""""""""""""""""""""""""""""""""
Groenheide et al. 2026
Barrovian metamorphism signals syn-orogenic delamination


HF_q0BC:
Heat Flux q = 0 boundary condition
No incoming heat flux and temperature varies freely at the lower boundary 
"""""""""""""""""""""""""""""""""""""""""""""""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

#%%
###########################################################
# Experiment selection and control parameter input.
###########################################################

#Experiment 1: Isothermal Relaxation (IR)
#Experiment 2: Radioactive Heating (Rad)
#Experiment 3: Mantle Heating  (MH)
#Experiment 4: Lithospheric Unzipping (LU)

experiment = 2
print('experiment=',experiment)

#Selecting the age and depth for the plots
plot_age = 30      # Myr
plot_depth = 30   # km (Interface crystalline basement-sedimentary cover = 30 km)
print('Temperature-Time plot at z =', plot_depth, 'km')
print("Temperature-Height plot at t =", plot_age, 'Myr')


# grid resolution
Nz = 101

# layer thicknesses (m)
L1 = 25e3
L2 = 5e3
L3 = 25e3

# temperature ranges (K)
T_moho = np.linspace(400,900,6) + 273.15
T_mantle = np.linspace(1100,1400,4) + 273.15
T_surface = 273.15



#%%
###########################################################
# Creating the folder structure.
###########################################################

if experiment == 1:
    data_folder = "IR_data"
    fig_folder  = "IR_fig"

elif experiment == 2:
    data_folder = "Rad_data"
    fig_folder  = "Rad_fig"

elif experiment == 3:
    data_folder = "MH_data"
    fig_folder  = "MH_fig"

elif experiment == 4:
    data_folder = "LU_data"
    fig_folder  = "LU_fig"

else:
    sys.exit("Experiment not valid")

# create folders if they do not exist
os.makedirs(data_folder, exist_ok=True)
os.makedirs(fig_folder, exist_ok=True)

#%%
###########################################################
#Building the grid.
###########################################################

# Layer 1: Crystalline basement (Diorite)
# Layer 2: Sedimentary cover (Carbonates)
# Layer 3: Overburden layer (Carbonates)
# ---
#  | L3
#  | L2
#  | L1
# ---


def build_grid(L1,L2,L3,Nz):
    L=L1+L2+L3                      #Height of domain                    
    z = np.linspace(0,L,Nz)
    h = L/(Nz-1)                    #Cell size
    return z, h

#%%
###########################################################
#Material Parameters.
###########################################################

# rho: density (kg/m3)
# hcond: heat conductivity (W/mK)
# hcapa: heat capacity (J/kgK)
# kappa: Thermal diffusivity of layer (m2/s)
# hprod: Radioactive heat production (µW/m3)

def material_properties(z,L1,L2):

    Nz = len(z)

    rho = np.zeros(Nz)
    hcond = np.zeros(Nz)
    hcapa = np.zeros(Nz)
    hprod = np.zeros(Nz)
    kappa = np.zeros(Nz)

    if experiment in [2,4]:             #Implementing radioactive heating in experiment 2 and 4.
        hprod1 = 1e-6
    else:
        hprod1 = 0
        
    for i in range(Nz):

        if z[i] < L1:                   #Crystalline basement
            rho[i] = 2870
            hcond[i] = 3.27
            hcapa[i] = 1140
            hprod[i] = hprod1

        elif z[i] < L1+L2:              #Sedimentary cover
            rho[i] = 2765
            hcond[i] = 2.12
            hcapa[i] = 780

        else:                           #Overburden
            rho[i] = 2765
            hcond[i] = 2.12
            hcapa[i] = 780

        kappa[i] = hcond[i] / (rho[i]*hcapa[i])
    
    return rho, hcond, hcapa, hprod, kappa

#%%
###########################################################
#Time-related parameters.
###########################################################
def time_parameters(experiment, h, kappa):

    year = 365.25 * 24 * 3600
    Myr  = 1e6 * year
    CFLnb = 0.25                            #Courant-Friedrichs-Lewy number 
    kappa_max = np.max(kappa)
    
    dt = h**2/(2*kappa_max) * CFLnb
    
    
    #Selecting the runtime for the experiments
    if experiment in [1,2]:
        tfinal = 30e6 * year                
        nstep  = int(tfinal/dt)
        if plot_age > tfinal/Myr:
            raise ValueError("Plot age outside simulation runtime")
    elif experiment in [3,4]:
        tfinal = 10e6 * year
        nstep  = int(tfinal/dt)
        if plot_age > tfinal/Myr:
            raise ValueError("Plot age outside simulation runtime")
    else:
        raise ValueError("Invalid experiment")

    return year, Myr, dt, tfinal, nstep

#%%
###########################################################
#Constructing the initial temperature profiles.
###########################################################

def initial_temperature(z,L1,L2,L3,T_surface,T_moho):

    T = np.zeros(len(z))

    for i in range(len(z)):

        if z[i] < L1+L2:
            T[i] = z[i]/(L1+L2)*(T_surface-T_moho) + T_moho

        else:
            T[i] = (z[i]-(L1+L2))/L3*(T_surface-T_moho) + T_moho

    return T

#%%
###########################################################
#Boundary conditions
###########################################################
def boundary_conditions(Tmantle, T_surface):

    BC_top = T_surface

    if experiment in [1, 2]:
        #No fixed temperature and no incoming heat flux.
        BC_bttm = 0

    elif experiment in [3, 4]:
        # Fixed temperature at base.
        BC_bttm = Tmantle

    return BC_top, BC_bttm


#%%
###########################################################
#The actual thermal model.
###########################################################

def run_model(T_initial,z,h,dt,nstep,
              rho,hcond,hcapa,hprod,
              BC_top,BC_bttm):

    Nz = len(z)

    T_old = T_initial.copy()
    T_new = np.zeros(Nz)

    results = []

    time = 0

    for istep in range(nstep):

        for i in range(1,Nz-1):

            T_new[i] = T_old[i] \
                + dt/(rho[i]*hcapa[i]*h**2) * (
                    (hcond[i]+hcond[i+1])/2*(T_old[i+1]-T_old[i])
                    - (hcond[i]+hcond[i-1])/2*(T_old[i]-T_old[i-1])
                ) \
                + hprod[i]/(rho[i]*hcapa[i])*dt
                
                
        if experiment in [1, 2]:
            #Completely skip boundary condition and let it cool freely.
            T_new[0] = T_new[1] 
        elif experiment in [3, 4]:
            T_new[0] = BC_bttm
        
        #fix surface T
        T_new[-1] = BC_top
        
        T_old[:] = T_new[:]

        time += dt

        results.append(T_new.copy())

    times = np.arange(1,nstep+1)*dt    
    return np.array(results), times

#%%
###########################################################
#Saving our results.
###########################################################
def save_results(folder,filename,time,results,z):

    path = os.path.join(folder,filename)

    with open(path,"w") as f:

        header = "Time(Ma)\t" + "\t".join([f"{d/1000:.2f}" for d in z])
        f.write(header+"\n")

        for i in range(len(time)):
            row = [time[i]] + list(results[i])
            f.write("\t".join(map(str,row))+"\n")

#%%
###########################################################
#Data visualisation
###########################################################
#Temperature-depth plot.

def plt_Tz(data, z, age, filename):
    folder = fig_folder
    
    plt.figure(dpi=250)
    
    #Looping over all runs eg. initial conditions.
    for run in data:                       

        results = run["results"]
        time = run["time"]
        label = run["label"]

        t_index = np.argmin(np.abs(time - age))

        T = results[t_index] 

        plt.plot(z/1000, np.flip(T-273.15), label= label)

    plt.xlabel("Depth (km)")
    plt.ylabel("Temperature (°C)")
    plt.gca().invert_xaxis()
    plt.axhline(y=700, color='black', linestyle='--', linewidth=0.7)
    plt.axhline(y=550, color='red', linestyle='--', linewidth=0.7)
    plt.axvline(x=plot_depth, color='black', linestyle='--', linewidth=0.7)
    plt.grid(True)
    plt.legend()
    plt.title(f"Temperature against Depth at {age/(1e6*365.25*24*3600):.1f} Myr")
    plt.tight_layout()
    plt.savefig(os.path.join(folder,filename))
    plt.show()
    
#Temperature-time plot.   
def plt_Tt(data, z, depth, filename):
    
    folder = fig_folder
    
    plt.figure(dpi=250)

    depth_index = np.argmin(np.abs((L1+L2+L3 - z) - depth*1000))

    for run in data:                        

        results = run["results"]
        time = run["time"]
        label = run["label"]

        T_time = results[:, depth_index] 

        plt.plot(time/(1e6*365.25*24*3600), T_time-273.15, label=label)

    plt.xlabel("Time (Myr)")
    plt.ylabel("Temperature (°C)")
    plt.legend()
    plt.axhline(y=700, color='black', linestyle='--', linewidth=0.7)
    plt.axhline(y=550, color='red', linestyle='--', linewidth=0.7)
    plt.grid(True)
    plt.title(f"Temperature Evolution at {depth} km Depth")
    plt.tight_layout()
    plt.savefig(os.path.join(folder,filename))
    plt.show()
    
#%%
###########################################################
#Combining the modules into a main program.
###########################################################

def main():
    
    #Building the grid.
    z, h = build_grid(L1,L2,L3,Nz)                          

    #Defining the material properties.
    rho,hcond,hcapa,hprod,kappa = material_properties(z,L1,L2)
    
    #Defining the time_parameters.
    year,Myr,dt,tfinal,nstep = time_parameters(experiment,h,kappa)

    ###########################################################
    #Experiment 1 and 2
    ###########################################################
    if experiment in [1,2]:
        data_storage = []
        for Tmo in T_moho:

            print(f"Calculating M_{round(Tmo-273.15,1)}")

            #Constructing initial temperature profiles.
            T0 = initial_temperature(z,L1,L2,L3,T_surface,Tmo)

            #Boundary conditions: fixed flux at base, fixed T at surface.
            BC_top, BC_bttm = boundary_conditions(T_mantle[0], T_surface)

            #Running the actual model.
            results, time = run_model(
                T0,z,h,dt,nstep,
                rho,hcond,hcapa,hprod,
                BC_top, BC_bttm
            )
            
            #Saving the results.
            save_results(data_folder,f"M_{round(Tmo-273.15,1)}.ascii",time,results,z)
            
            #Locally storing data for visualisation.
            data_storage.append({"results": results,
                "time": time,
                "label": f"{round(Tmo-273.15,1)}°C"
            })
           
        # Temperature at fixed time (depth plot)
        fn_Tz=f'Tz_{plot_age}.png'
        plt_Tz(data_storage, z, plot_age * Myr, fn_Tz)
        
        # Temperature at fixed depth (time plot)
        fn_Tt=f'Tt_{plot_depth}.png'
        plt_Tt(data_storage, z, plot_depth, fn_Tt)
    
    
    ###########################################################
    #Experiment 3 and 4
    ###########################################################
    elif experiment in [3,4]:
        data_FMo = []
        data_FMa = []
        
        for Tmo in T_moho:
            T0 = initial_temperature(z,L1,L2,L3,T_surface,Tmo)
            for Tma in T_mantle:
                
                print(f"Calculating M_{round(Tmo-273.15,1)}_{round(Tma-273.15)}")
                
                BC_top, BC_bttm = boundary_conditions(Tma, T_surface)
                results, time = run_model(
                    T0,z,h,dt,nstep,
                    rho,hcond,hcapa,hprod,
                    BC_top, BC_bttm
                )
                
                save_results(data_folder,f"M_{round(Tmo-273.15,1)}_{round(Tma-273.15)}.ascii",time,results,z)
                
                if np.isclose(Tmo-273.15, 900):
                    data_FMo.append({"results": results,
                        "time": time,
                        "label": f"{round(Tma-273.15,1)}°C"
                    })
                
                if np.isclose(Tma-273.15, 1400):
                    data_FMa.append({"results": results,
                        "time": time,
                        "label": f"{round(Tmo-273.15,1)}°C"
                    })
                    
        fn_Tz_Ma=f'Tz_{plot_age}_Ma{round(Tma-273.15,1)}.png'
        fn_Tz_Mo=f'Tz_{plot_age}_Mo{round(Tmo-273.15,1)}.png'
        plt_Tz(data_FMo, z, plot_age * Myr, fn_Tz_Ma)
        plt_Tz(data_FMa, z, plot_age * Myr, fn_Tz_Mo)
        
        fn_Tt_Ma=f'Tt_{plot_depth}_Ma{round(Tma-273.15,1)}.png'
        fn_Tt_Mo=f'Tt_{plot_depth}_Mo{round(Tmo-273.15,1)}.png'
        plt_Tt(data_FMo, z, plot_depth, fn_Tt_Ma)
        plt_Tt(data_FMa, z, plot_depth, fn_Tt_Mo)

###########################################################
#Running everything.
###########################################################
if __name__ == "__main__":
    main()