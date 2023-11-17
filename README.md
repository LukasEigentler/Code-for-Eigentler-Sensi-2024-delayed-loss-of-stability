# Code-repository-delayed-loss-of-stability

This repository contains Matlab and AUTO code that was used to produce the results shown in the paper L. Eigentler, M. Sensi (2023). Delayed loss of stability of periodic travelling waves: insights from the analysis of essential spectra.

All code is provided through GNU General Public License v3.0. Please see the license file for more information

The folder contains two folders: "Klausmeier" and "Mussels". Each folder contains code for the respective model presented in the paper. Below, I outline the contents for the "Klausmeier" folder. The contents of the "Mussels" folder is identical (different filenames).

The "Klausmeier" folder contains two subfolders: "Num_sim" and "Num_cont". 

In the "Num_sim" folder, there are the following files:
- "klausmeier_single_sim.m" solves the system once.
- "klausmeier_delay_pred.m" performs a prediction of the delay and numerically solves the system to make a comparison.
- "klausmeier_memorylessness.m" performs a range of simulations to visualise the memorylessness property of the delayed loss of stability phenomneon with respect to the dynamics that occur before crossing the stability boundary.
- "klausmeier_time_delay_rate_of_change.m" solves the system for a parameter change regime in which the bifurcation parameter decreases at a linear rate. It performs simulations for a range of different decay rates, analyses the data and loads data from corresponding essential spectra for comparison between solution data and PTW spectra.
- "klausmeier_time_delay_sudden_change_A.m" solves the system for a parameter change regime in which the bifurcation parameter decreases instantaneously to a target value. It performs simulations for a range of different dtarget values, analyses the data and loads data from corresponding essential spectra for comparison between solution data and PTW spectra.
- "importsol_klausmeier.m" is a Matlab function to import solution data held in "s." files from AUTO.
- "importsolutiondata_klausmeier.m" is a Matlab function to import data held in "b." files from AUTO.
- "klausmeierode.m" is a Matlab function that contains the equation data for use in Matlab's ode solver.

The "Num_cont" folder contains the following files:
- "stab_diag_klausmeier.m" loads numerical continuation data and visualises the Busse balloon.
- "parameters.dat" stores the model's parameter value (see accompanying readme file) for use in numerical continuation.

Futher, the "Num_cont" folder contains subfolders with the following files:
"Pattern_generation" folder:
- "pattern_gen.auto" generates a PTW through numerical continuation.
- "wavelength_cont.auto" calculates wavelength contours for the Busse balloon.
- "fold_klausmeier.auto" calculates the location of a fold in the Busse balloon.
- "homoclinic_klausmeier.auto" calculates the location of a homoclinic orbit (approx by large wavelength) in the Busse balloon.
- "hopf_klausmeier.auto" calculates the location of a Hopf bifurcation in the Busse balloon.
- "importhomoclinic_klausmeier.m" is a Matlab function to import the data on homoclinic orbits into Matlab.
- "importwavelength_klausmeier.m" is a Matlab function to import the data on wavelength contours into Matlab.

"Pattern_stability" folder:
- "spectrum_cont_klausmeier.auto" uses numerical continuation to calculate an essential spectrum from externally supplied starting data.
- "gamma_0_klausmeier.m" solves the Matrix eigenvalue problem that is subsequentally fed into the numerical continuation to calculate an essential spectrum or stability boundaries.
- "gamma_0_klausmeier_spectra_batch_calc.m" is as above but automated to be called from a bash file as part of the spectrum calculation.
- "mult_spectra_plot.m" visualises data from multiple spectra.
- "spectrum_post_proc_batch.m" analyses spectra.
- "import....m" files import data as suggested by the filename into Matlab.
- "Eckhaus_stab_boundary/eckhaus.auto" detects and continues Eckhaus stability boundaries.

"Spectra_batch_calc" folder:
- "spectra_batch_bash.bash" is a bash script that calculates multiple spectra.

 
