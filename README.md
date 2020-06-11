#  Using the LIST model to Estimate the Effects of Contact Tracing on COVID-19 Endemic Equilibria in England and its Regions



By Rosalyn J. Moran, Alexander J. Billig, Maell Cullen, Adeel Razi, Jean Daunizeau, Rob Leech, Karl J. Friston


Governments across Europe are preparing for the emergence from lockdown, in phases, to prevent a resurgence in cases of COVID-19. Along with social distancing (SD) measures, contact tracing – find, track, trace and isolate (FTTI) policies are also being implemented. Here, we investigate FTTI policies in terms of their impact on the endemic equilibrium. We used a generative model – the dynamic causal ‘Location’, ‘Infection’, ‘Symptom’ and ‘Testing’ (LIST) model to identify testing, tracing, and quarantine requirements. We optimised LIST model parameters based on time series of daily reported cases and deaths of COVID-19 in England—and based upon reported cases in the nine regions of England and in all 150 upper tier local authorities. Using these optimised parameters, we forecasted infection rates and the impact of FTTI for each area—national, regional, and local. Predicting data from early June 2020, we find that under conditions of medium-term immunity, a ‘40%’ FTTI policy (or greater), could reach a distinct endemic equilibrium that produces a significantly lower death rate and a decrease in ICU occupancy. Considering regions of England in isolation, some regions could substantially reduce death rates with 20% efficacy. We characterise the accompanying endemic equilibria in terms of dynamical stability. These analyses suggest that FTTI will not only save lives, it could underwrite the stability of any endemic steady-state we manage to attain.



Code used to generate the results and figures in the paper in scripts:

SPM12 is required to reproduce these results (https://www.fil.ion.ucl.ac.uk/spm/covid-19/#software). 

Version SPM12 Toolbox files used for this work found in spm_version_files

To Run Model:

Add spm to matlab path

Then

To Run England: DEM_COVID_England.m
To Run Region:  DEM_COVID_NorthWest.m

To study Bifurcation: Set Ep_Nation.ttt      

For Example to Test Endemic Effects of an 80% successful Find Track and Trace policy

Ep_Nation =  log(0.80); 

