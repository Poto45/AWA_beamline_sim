# AWA_beamline_sim
This repository will be for the full beamline simulation that was started in Fall 2023 for the experiment in 2024.
The method of doing a full beamline simulation is:
1. Opal-T for gun to Yag4, because nothing in that section changes; optimized to minimize $\epsilon_x$, allowing the beam to pass through the small structure, while trying to reduce $s_{rms}$, creating a higher peak current. Increasing peak current (desired for creating wakefields) has the potential to blow up emittance, hence the inclusion of both. 
2. Elegant simulation from Yag4 to the structure to optimize quadrupoles to create beam waist in center of structure
3. Take optimized quadrupole information and put into Opal-T for Yag4 to start of structure
4. WarpX to propagate beam (uploaded from Opal-T) through the structure
5. Elegant to optimize quadrupoles from end of structure to create beam waist at final Yag screen 
6. Opal-T from end of structure to final Yag screen using optimized quadrupole settings

** Sirepo was used for the Elegant portions and the information is also uploaded here. **

** Opal-T and WarpX were run on Bebop/Swing, respectively; queue scripts are also uploaded here. **

** Also created, a toy model of the Longitudinal Phase Space (LPS) measurement system. This is to just better understand the modulations and distortions of LPS. **

## Opal-T


## Elegant
Input is SDDS file, but you also need the numerical value of the Central momentum of the beamline [MeV/c]. For individual particles' information, x, y, z are in [m], compared to where the reference particle are and the momenta are in $\beta \gamma$. These are the same units as used in Opal-T. 
Process for Opal-T to Elegant:
1. Opal-T outputs an h5 file with all information. 

## WarpX





