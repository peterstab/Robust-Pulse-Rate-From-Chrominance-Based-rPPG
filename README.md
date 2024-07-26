# Robust-Pulse-Rate-From-Chrominance-Based-rPPG
Matlab code for the algorithm described in the paper Robust_Pulse_Rate_From_Chrominance_Based_rPPG(1).pdf


%INDICATIONS

%If you want to consult the RGB and reference data:

% 'starting_workspace.mat': this workspace contains the starting RGB and
% reference data. 
%However, they are private data. For each subject we worked with (67850 x 3) double matrix,
% meaning approximately 10 minutes (sampling frequency 115 Hz) of R, G & B signals. Reference data are (67850 x 1) double.

%The only scripts to be run are the following two:

%'general_procedure.m': through this script you can follow each passage
%that we applied in order to work on the starting data
%and finally extract the pulse rate. 
%At the beginning, you have to select a single patient in the variable
%'patient' and the phase in the variable 'segment' that you want to analyse, and
%the results will be provided for this specific choice. In particular, the
%plots of the signals, the plots of the spectra, the SNR and the BPM are 
%provided for each method.

%'global_results.m': the workspace will be loaded and some indexes (RMSE, STD, 
%Pearson, slope of linear fit) are computed for each method and for each phase
%over the 20 patients, using a function called 'computing_data.m'. This function has
%the same code seen in 'general_procedure.m', but with the necessary 
%modifications to provide all data for 20 patients without every single plot.

%The other functions present in the folder are automatically called by 
%'general_procedure.m' and 'global_results.m'.
