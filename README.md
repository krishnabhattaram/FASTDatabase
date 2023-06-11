# FAST Database
> [Database Location][Database]

:warning: This README is currently being updated as the database is in flux. Please take careful note of the documentation to ensure the data product is as expected, and reach out to the project contacts listed above with any doubts

### Author Information:
- Principal Investigator: [Chris Chaston][ccc], [Krishnakumar Bhattaram][kkb], Ennio Sanchez
- Points of Contact: [Chris Chaston][ccc], [Krishnakumar Bhattaram][kkb]



## Introduction
The Fast Auroral Snapshot Explorer (FAST) spacecraft collected up to five gigabytes per day over its operational lifespan, consolidating a rich variety of particle and fields data over a range of events and altitudes. Among other uses, these data have proven invaluable for application to the study and modeling of wave-particle interactions in Alfvénic and other waves modes, ionospheric energization, and ion outflow. However, instrumentation and operational idiosyncrasies and outdated access and representation protocols have proven to be a barrier to their usage for researchers outside the FAST team and associates. We have therefore compiled a consolidated, sanitized, and standardized database of FAST particle, field, and spectral energy quantities on a common time axis to facilitate easy access for pursuing these scientific goals without requiring detailed knowledge of instrument operation modes and instrument characteristics. Quantities in the database include spin-period averaged vector electric and magnetic field magnitudes and spectral energy densities from the sub-Hz to kHz range, electron and ion number and energy in- and outflows, and geospatial positioning and trajectory data, which are compiled as ISTP/IACG-compliant Common Data Format (CDF) files. These data are made publicly available [here][Database]. Representations of the reduced measurements comprising the database are presented from the beginning of the mission in 1996 to the failure of the electric field instrument in the early 2000s. 

## Database Overview

<p align="center">
 <img width="529" alt="image" src="https://github.com/krishnabhattaram/FASTDatabase/assets/60773678/4d1567f7-5a50-4f87-b9b4-b1aa4a35375c">
</p>

The data collected in the FAST Database is sourced from the Satellite Data Tool ([SDT]) developed at the UC Space Sciences Laboratory. The data mutation process is highlighted in the above diagram. The DC vector despun magnetic field perturbation is extracted through the UCLA despin procedure, while despun DC electric field quantities are derived from SSL's fields despin procedure using the V5-V8 and V1-V2 probes ([Despin Sources]). Spectral energy densities are also sourced from the FAST onboard Digital Signal Processing (DSP) instrument. Particle data are presented integrated over all or a quartile of pitch angles, of which the exact procedure is described in the below particles section. Fields and particle quantities are averaged over the FAST five-second spin interval. Telemetry data is presented at the same cadence. Further details for the database quantities are avilable under the 'catdesc' attribute of the 'vatt_general' field in the CDF metadata per quantity.

### Telemetry Quantities
Attribute Name | Units |	Source Instrumentation |	Description 
--- | --- | --- | --- 
b_foot/b_model	| nT |	Telemetry/IGRF |	Footpoint magnetic field (3-component), Magnetic field from model (3-component) in Cartesian Coordinates
fa_pos |	km |	Telemetry | 	Spacecraft position (3-component, GEI)
lat	| deg	| Telemetry |	Spacecraft latitude
lng	| deg	| Telemetry	| Spacecraft longitude
alt	| km |	Telemetry |	Spacecraft altitude
flat | deg | Telemetry |	Magnetic footpoint latitude
flng |	deg	 | Telemetry |	Magnetic footpoint longitude
ilat |	deg	| Telemetry |	Invariant latitude
ilng |	deg |	Telemetry |	Invariant longitude
mlt |	hours |	Telemetry |	Magnetic local time
fa_vel |	km/s |	Telemetry |	Spacecraft velocity (3-component, GEI)

Quantities for telemetry are sourced from the tplot quantities made available from the despin procedures. 
### Fields Quantities

<p align="center">
 <img width="510" alt="image" src="https://github.com/krishnabhattaram/FASTDatabase/assets/60773678/1b136bb7-69dd-4994-bd38-0fee129fa53b">
</p>

The magnetic field measurements are sourced from three separate data quantities from the FAST satellite; the fluxgate magnetometer, search coil magnetometer, and the on-board DSP instrument. The spectral energy quantities are calculated through a fast-fourier transform on the fluxgate and seach-coil data (magz-sp and mag3ac-sp respectively), with the DSP quantity likewise providing the results of an onboard FFT. Note that the search coil magnetometer spectral range is 0-50 Hz, while the search coil is rated between 50 Hz to 16 kHz in its single-axis mode. Users are encouraged to validate the spectral range of the data (and account for noise floor limitations of the search coil data). Further details are available at [FAST Fields][Fields Source].

The electric field instrument radial field booms carrying probes 3 and 4 did not deploy or did not deploy fully; for this reason, the data presented in this database use measurements from probes 1/2 and 5/8. The electric field instrument also developed an intermittent DC offset after January 2000. This database therefore does not provide data beyond this date. The electric field instrument has also been observed to generate data unphysical artifacts as the instrument turns on approaching the auroral oval. In an attempt to alleviate this, *all electric field measurements beyond 1 V/m are discarded*. Please take these modifications into account, and contact the contacts of the project if any questions arise. 

Attribute Name | Units |	Source Instrumentation |	Description |	Frequencies 
--- | --- | --- | --- | --- 
magz-sp	| nT<sup>2</sup>/Hz	| Fluxgate Magnetometer	| Z-axis s/c spin axis spectral energy density from fluxgate magnetometer	| 1 to 160 Hz, steps of 1.25 Hz
mag3ac-sp	| nT<sup>2</sup>/Hz	| Search Coil Magnetometer	| Instrument (Mag3ac) z-axis search coil spectral energy density	| 1 to 160 Hz, steps of 1.25 Hz
eav-sp |(mV/m)<sup>2</sup>/Hz| Electric Field Instrument |	Electric field spectral energy density along spacecraft trajectory (capped at 1 V/m) |	1 to 160 Hz, steps of 1.25 Hz, plus [0.05, 0.1, 0.2, 0.5]
magz/mag3ac/eav-sp-tres |	s	| FGM/Electric Field/SCM | Time resolution of FGM z-axis/SCM z-axis/electric field time-series measurements prior to entering the spectral transform logic | -
dspV5v8-sp | log((V/m)<sup>2</sup>/Hz)	| Electric Field Instrument |	Spin plane electric field spectral density. Calculated through FFT through spacecraft DSP (capped at 1 V/m)	| Log-scale, 16 frequencies from 3.16 to 17800
dspmag3ac-sp| log(nT<sup>2</sup>/Hz) | Search Coil Magnetometer | Search coil z-axis spectral energy density. Calculated through FFT through spacecraft DSP | 	Log-scale, 16 frequencies from 3.16 to 17800
e_avg |	mV/m | Electric Field Instrument 	| Spin averaged electric field along spacecraft trajectory | -
e_sdev | mV/m |	Electric Field Instrument |	Standard deviation in total electric field per spin interval	| -
db_avg |	nT |	SCM/FGM	| dB (rel. to IGRF) along trajectory, dB cross track, dB field-aligned |	-
tot_b_avg	| nT	| SCM/FGM |	Spin Period Averaged Total Magnetic Field |	-
tot_b_sdev |	nT | SCM/FGM |	Standard deviation in total magnetic field per spin interval |	-

### Particles Quantities
Two variants of measurements are provided in this database; the first includes number/energy flux of the ions and electrons integrated over the entire 2π radians of pitch angle, while the second only integrate over a single quandrant. The former quantity, titled "outflow" or "inflow" quantities, provide the ion number and energy flux in the out- or inflow direction opposite to the spacecraft motion (to avoid ram effects). Energy ranges were selected with a focus on studying outflow effects. 

Attribute Name |	Units	| Source Instrument |	Description |	Integration Angles |	Energies
--- | --- | --- | --- | --- | --- 
j_2d_fa_ees |	1/s-cm2 |	ESA	| Electron number flux 	| All Pitch Angles |	20 eV to 30 keV
je_2d_fa_ees |	mW/s-m2 |	ESA	| Electron energy flux | 	All Pitch Angles |	20 eV to 30 keV
j_2d_fa_ies |	1/s-cm2	| ESA |	Ion number flux | All Pitch Angles |	20 eV to 1 keV
je_2d_fa_ees | mW/s-m2 |	ESA	| Ion energy flux |	All Pitch Angles |	20 eV to 1 keV
j_2d_fa_ies_outflow and inflow | 1/s-cm2 | ESA | Ion number flux |	In/Outflow quadrant opposite sc_vel | 	20 eV to 1 keV
je_2d_fa_ies_outflow and inflow	| mW/s-cm2 | ESA | Ion energy flux | In/Outflow quadrant opposite sc_vel | 20 eV to 1 keV

## Code Description
All code was run with SDT configured on the machine as an SDT batch job. lws_stats_sanitized.pro is the code that generates the main data products present in the database. lws_stats_supplement_data.pro includes time resolution information. Both files make use of the helper functions located in lws_stats_helper_functions.pro. Refer any questions to the contacts of the project. 

   [ccc]: <mailto: cchaston@berkeley.edu>
   [kkb]: <mailto: krishnabhattaram@berkeley.edu>
   [Database]: <https://drive.google.com/drive/u/1/folders/1mLdcruLAfEBWqEi9wCFRof3kqzO-VA3i>
   [SDT]: <http://sprg.ssl.berkeley.edu/fast/scienceops/fastidl.html>
   [Despin Sources]: <http://sprg.ssl.berkeley.edu/fast/scienceops/fast_fields_help.html>
   [Fields Source]: <http://sprg.ssl.berkeley.edu/fast/scienceops/fast_fields_help.html>
