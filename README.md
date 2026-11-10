**Groenheide et al. 2026**

**Barrovian metamorphism signals syn-orogenic delamination.** 

**Read me:**
Exp. 1 & 2:
CTBC: Fixed Moho temperature at the base of the system 
HF_q0BC: Letting temperature vary freely and no incoming heat flux at the base. 
HF_qcBC: Letting temperature vary freely and having a constant heat flux of 25 mW/m2.
HF_LAB_Rob: Letting temperature vary freely and having a constant heat flux based on the temperature difference between Moho and Lithosphere-Asthenosphere Boundary. 

Exp 3 & 4:
HF_LAB_Rob: Using Robin boundary conditions as a contact between the crust (Moho) and a hot astehonsphere directly below. Formation of new mantle lithosphere is modelled as a growing half-sapce thermal boundary layer. 

Plotting:
Exp1_2_BC_comparison: This code uses the data of all boundary conditions to plot them for a given temperature allowing for comparison (in our case at a Moho T of 700 ◦C ).

HF_T_t: This code compares different heating windows for a fixed moho-mantle temperature combination (in our case 700-1300 ◦C ). It plots both the Heat flux (calculated from the temperatures) and temperature against time.

Conclusion_plot: This plot uses our data from exp 2 and 4 to be able to compare processes. It plots the end-member temperature-time profiles at our location of interest. 

**Background:**
The codes used for this paper consist in two categories: Modelling and Plotting. All Modelling codes use the same structure. For the different Boundary Conditions in experiment 1 and 2 we have different codes, varying only in BC's. Exp 3 and 4 boundary conditions were still a fixed temperature at the base as a remnant of a first starting code. For Experiment 3 and 4 we used HF_LAB_Rob where we uses a Robin boundary condition together with a growing half-spacce thermal boundary condition to model direct asthenossphere contact with a new growing mantle litosphere. Visualisation was either done within the model code or with separate codes. 
