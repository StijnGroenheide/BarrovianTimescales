# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 20:29:19 2025

@author: Stijn
"""

"Guided Research: Data processing conclusion plot"

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

#%% 
###########################################################
# Data selection: Experiment 2: radioactive heating model
###########################################################



# Path to the folder containing 
txt_folder_rad = "Rad_data"

# Getting the names of the files
file_names_rad = [f for f in os.listdir(txt_folder_rad) ]

 
for file_name in file_names_rad:
    df_name = file_name[:-6] # Remove ".ascii"
    globals()[df_name] = pd.read_csv(f"{txt_folder_rad}/{file_name}", sep='\t', header=0)


datasets_rad = []


for file_name in file_names_rad:
    df_name = file_name[:-6]   
    df = pd.read_csv(f"{txt_folder_rad}/{file_name}", sep='\t', header=0)
    
    # Store as a variable
    globals()[df_name] = df  
    
    datasets_rad.append((df_name, df))  

#%%
###########################################################
# Data selection: Experiment 4: Combined heating model
###########################################################



txt_folder = "LU_data_30Myr"

file_names_LU = [f for f in os.listdir(txt_folder)]

for file_name in file_names_LU:
    df_name = file_name[:-6]   
    globals()[df_name] = pd.read_csv(f"{txt_folder}/{file_name}", sep='\t', header=0)

datasets_LU = []

for file_name in file_names_LU:
    df_name = file_name[:-6]  
    df = pd.read_csv(f"{txt_folder}/{file_name}", sep='\t', header=0)
    
    # Store as a variable
    globals()[df_name] = df   
    datasets_LU.append((df_name, df))  

#%% 
###########################################################
#Processing and seperating the data in groups for comparison
###########################################################


#Creating groups based on Moho T
Mo_400 = []
Mo_500 = []
Mo_600 = []
Mo_700 = []
Mo_800 = []
Mo_900 = []



Mo_group = [Mo_400, Mo_500, Mo_600, Mo_700, Mo_800, Mo_900]

T_labels_mo = []
for df_name, df in datasets_LU:
    parts = df_name.split("_")  # Splitting it up in parts
    moho_value = parts[1]  # Moho value
    T_labels_mo.append(parts[1]) 
    
    if moho_value == "400.0":
        Mo_400.append((df_name, df))
    elif moho_value == "500.0":
        Mo_500.append((df_name, df))
    elif moho_value == "600.0":
        Mo_600.append((df_name, df))
    elif moho_value == "700.0":
        Mo_700.append((df_name, df))
    elif moho_value == "800.0":
        Mo_800.append((df_name, df))
    elif moho_value == "900.0":
        Mo_900.append((df_name, df))

#Creating groups based on Mantle temperature 
Ma_1100 = []
Ma_1200 = []
Ma_1300 = []
Ma_1400 = []

Ma_group = [Ma_1100, Ma_1200, Ma_1300, Ma_1400] 


T_labels_ma = []
for df_name, df in datasets_LU:
    parts = df_name.split("_")  
    ma_value = parts[2]     #extracting mantle temperatures
    T_labels_ma.append(parts[2]) 
    if ma_value == "1100":
        Ma_1100.append((df_name, df))
    elif ma_value == "1200":
        Ma_1200.append((df_name, df))
    elif ma_value == "1300":
        Ma_1300.append((df_name, df))
    elif ma_value == "1400":
        Ma_1400.append((df_name, df))
        
#print(Ma_1100)


#%%
###########################################################
#Plotting
###########################################################

fixed_height_km = 25 #= 30 km depth
H = np.linspace(0, 55, 101)
height_idx = (np.abs(H - fixed_height_km)).argmin()

# Time axis in Myr
time_ma = np.linspace(0, 30, np.size(df,0))

plt.figure(dpi=250)

###########################################################
#Radioactive heating model plot
###########################################################

for df_name, df in [datasets_rad[0]]: #Selecting the minimum initial conditions: 400 
    parts = df_name.split('_')
    moho_value = int(float(parts[1]))

    T_series = df.iloc[:, height_idx + 1]  # +1 because first column is time
    plt.plot(time_ma, T_series - 273.15, color = 'blue', linestyle = 'dashed', label = f'{moho_value}°C')

for df_name, df in [datasets_rad[3]]: #Selecting the middle: 700 
    parts = df_name.split('_')
    moho_value = int(float(parts[1]))

    T_series = df.iloc[:, height_idx + 1]
    plt.plot(time_ma, T_series - 273.15, color = 'blue', label = f'{moho_value}°C')

for df_name, df in [datasets_rad[-1]]: #Selecting the maximum initial conditions:  900
    parts = df_name.split('_')
    moho_value = int(float(parts[1]))

    T_series = df.iloc[:, height_idx + 1]  # +1 because first column is time
    plt.plot(time_ma, T_series - 273.15, color = 'blue', linestyle = 'dashed', label = f'{moho_value}°C')

###########################################################
#Combined heating Model
###########################################################    
    
for df_name,df in [Ma_1100[0]]: #Here I now only select an upper mantle limit
    parts = df_name.split('_')
    moho_value = int(float(parts[1]))
    Ma_value = int(float(parts[2]))
    T_series = df.iloc[:, height_idx + 1]  # +1 because first column is time
    plt.plot(time_ma, T_series - 273.15, color='green', linestyle ='dashed', label = f'{moho_value}-{Ma_value}°C')
    
for df_name,df in [Ma_1300[3]]: #Here I now only select an upper mantle limit
    parts = df_name.split('_')
    moho_value = int(float(parts[1]))
    Ma_value = int(float(parts[2]))
    T_series = df.iloc[:, height_idx + 1]  # +1 because first column is time
    plt.plot(time_ma, T_series - 273.15, color='green', label = f'{moho_value}-{Ma_value}°C')
        
for df_name,df in [Ma_1400[5]]: #here I only select an lower mantle limit. 
    parts = df_name.split('_')
    moho_value = int(float(parts[1]))
    Ma_value = int(float(parts[2]))
    T_series = df.iloc[:, height_idx + 1]  # +1 because first column is time
    plt.plot(time_ma, T_series - 273.15, color='green', linestyle ='dashed', label = f'{moho_value}-{Ma_value}°C')  
    
    

    
# Plot formatting
plt.axhline(y=700, color='black', linestyle='--', linewidth=0.7)
plt.axhline(y=550, color='red', linestyle='--', linewidth=0.7)
plt.xlabel("Time (Myr)")
plt.ylabel("Temperature (°C)") 
plt.title("Temperature Range of the Radioactive and Combined Heating Model at 30 km")
plt.legend(title = 'Moho-Mantle Temperature', loc='upper center', bbox_to_anchor=(1.29,0.29), borderaxespad=0.0, ncol=2)
plt.grid(True)


