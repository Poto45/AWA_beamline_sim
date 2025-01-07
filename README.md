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
Currently, do not have a viable script to convert any distributions to Opal-T. 

Input units are:
[]
Output units are:
[$\beta_{x,y,z} \gamma$] 

## Elegant
Input is SDDS file, but you also need the numerical value of the Central momentum of the beamline [MeV/c]. For individual particles' information, x, y, z are in [m], compared to where the reference particle are and the momenta are in $\beta \gamma$. These are the same units as used in Opal-T. 
Input SDDS for "elegant" (as apposed to "spiffe") requires (x [m], xp [dimensionless], y [m], yp [dimensionless], t [s], p [$\beta \gamma$])
Process to Elegant:
1. Opal-T and WarpX outputs an h5 file with all information. This h5 file has a different heirarchy than what is allowed for the hdf2sdds used on the Argonne clusters, so we must change it to txt.
2. Use the 2Ele.py script to change it to a txt file. Changes are available for WarpX and Opal-T. This also prints the final momentum that is required to be manually input to Sirepo Elegant.
3. Finally, change the txt file to SDDS for Elegant, using: /lcrc/project/Bright-Beams/software/pelegant/elegantTree_osc/epics/extensions/bin/linux-x86_64/plaindata2sdds {txtfile} {newSDDSfilename} -inputMode=ascii "-separator= " -column=x,double,units='m' -col=xp,double -col=y,double,units='m' -col=yp,double -col=t,double,units='s' -col=p,double,units='m$be$nc'
4. To double check the SDDS conversion worked: /lcrc/project/Bright-Beams/software/pelegant/elegantTree_osc/epics/extensions/bin/linux-x86_64/sddsquery {SDDSfile}


## WarpX
I have Opal-T to WarpX file, but it needs to be commented out. No need to create Elegant to WarpX right now. Output units are: [$\gamma v$] for momenta. 




