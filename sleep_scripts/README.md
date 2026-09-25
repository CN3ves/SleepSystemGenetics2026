Main programs to extract sleep-wake/EEG phenotypes from *.smo files.

Files with the names of the mice, the group they belong to, and the number of mice per group are listed in 'Nrf1 directory.txt' and 'Drp1 Directory.txt' and read by the programs.

'FT' are the floxed/tamoxifen-injected mice (icKO), 'FV' the floxed/vehicle-injected mice (Tamoxifen controls), and 'CT' mice control/tamoxifen-injected mice (Genotype controls)

These Directory files also contain individual sleep onset after sleep deprivation (in 4s epochs afer light onset of day 3 with end of sleep deprivation = 5400) and 
the EEG reference used to normalize the EEG spectra to reduce variabilty in the inter-individual differences in EEG amplitude; in microVolt square).

So called *.smo files are binary files with following structure (a 'Record' in Pascal) for each 4 sec epoch in the 4-day experiment (86400 epochs total per smo file):
'''
 epoch = RECORD
            state          : char;
            bin            : array[0..400] of single;
            EEGv,EMGv,temp : single;
          END;
'''
'''
'state' can have the following values:
	'w','n','r' for EEG artefact-free waking, NREM sleep, and REM sleep, respectively.
	'1','2','3' for waking, NREM sleep, and REM sleep with EEG artefacts, respectively.
	'4','5','6' for waking, NREM sleep, and REM sleep with 'spindle-like' paroxysmal EEG activity, respectively.
	'9'         for Theta-Dominated Wakefulness (TDW).
'bin' the full spectal power densities at 0.25Hz resolution from 0-100Hz, i.e., 401 values, in volt-square per 0.25Hz 
'EEGv','EMGv','temp' for the EEG variance, EMG variance, and cortical temperature(set to 0.0 if not recorded).
'''