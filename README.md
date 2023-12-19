# AWA_beamline_sim
This repository will be for the full beamline simulation that was started in Fall 2023 for the experiment in 2024.
The method of doing a full beamline simulation is:
1. Opal-T for gun to Yag4, because nothing in that section changes
2. Elegant simulation from Yag4 to end of straight section to optimize quadrupoles to create beam waist in center of structure
3. Take optimized quadrupole information and put into Opal-T for Yag4 to start of structure
4. WarpX to propagate beam (uploaded from Opal-T) through the structure
5. Elegant to optimize quadrupoles from end of structure to create beam waist at final Yag screen 
6. Opal-T from end of structure to final Yag screen

*** Sirepo was used for the Elegant portions and the information is also uploaded here. ***
*** Opal-T and WarpX were run on Bebop/Swing, respectively; queue scripts are also uploaded here. ***

Also created, a toy model of the Longitudinal Phase Space (LPS) measurement system. This is to just better understand the modulations and distortions of LPS.







