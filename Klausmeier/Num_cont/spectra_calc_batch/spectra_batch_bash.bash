#!/bin/bash
# Bash script to automatically calculate spectra
# Author: Lukas Eigentler (lukas.eigentler@uni-bielefeld.de)
# License: GNU GPL
# Last updated: 17/11/2023
runind=0
while [ $runind -le 99 ]
do
	printf "Count has a value of $runind\n"
	((runind++))
	auto spectra_batch1.auto
	cd ../Pattern_stability/
	matlab -batch 'gamma_0_klausmeier_spectra_batch_calc'
	cd ../spectra_calc_batch/
	auto spectra_batch2.auto
	cd ../Pattern_stability/
	matlab -batch 'spectrum_post_proc_batch'
	cd ../spectra_calc_batch/
	runind++
done