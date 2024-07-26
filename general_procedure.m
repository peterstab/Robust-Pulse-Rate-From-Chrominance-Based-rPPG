%% load workspace

%This script was written in order to show each single processing stage
%for each algorithm.
clear all
load("starting_workspace.mat");

clearvars -except bvp total_data
close all




%% segment preparation

%Choose a patient between 1 and 20 using the variable 'patient':
patient=12;

%extraction of RGB data and reference data. We cut 20 seconds from data_tot
%due to an unbalance between reference data and RGB data. The 20 seconds
%were cut from the baseline phase.
data_tot=total_data{patient}(end-65549:end,:);
bvp_temp=bvp(2:end,patient);

%definition of sampling frequency
fs=115;

%Depending on the duration of each phase, a vector with indexes was defined
indexes=[fs*60*5 fs*60*7 fs*60*8 fs*60*8.5 fs*60*9 fs*60*9.5];

%Extraction of each phase from the total data
start_baseline=data_tot(1:indexes(1),:);
paced_breathing=data_tot(indexes(1)+1:indexes(2),:);
rest1=data_tot(indexes(2)+1:indexes(3),:);
hand_grip1=data_tot(indexes(3)+1:indexes(4),:);
rest2=data_tot(indexes(4)+1:indexes(5),:);
hand_grip2=data_tot(indexes(5)+1:indexes(6),:);
start_baseline_ref=bvp_temp(1:indexes(1),:);
paced_breathing_ref=bvp_temp(indexes(1)+1:indexes(2),:);
rest1_ref=bvp_temp(indexes(2)+1:indexes(3),:);
hand_grip1_ref=bvp_temp(indexes(3)+1:indexes(4),:);
rest2_ref=bvp_temp(indexes(4)+1:indexes(5),:);
hand_grip2_ref=bvp_temp(indexes(5)+1:indexes(6),:);
%The variable below was defined in order to extract only a portion of data
%from each phase, in order to find the portion with the smallest amount of
%interframe motion, according to the paper (2875/115= 25 seconds)
number_sample=2875;

%% segment choice
%Here you can select the phase changing the values of 'segment' with:
%start_baseline 
%paced_breathing 
%rest1
%hand_grip1 
%rest2
%hand_grip2

segment=start_baseline;

%'segment_ref' must correspond to the same phase chosen above (ex.
%start_baseline_ref if segment is 'start_baseline')

segment_ref=start_baseline_ref;

%% window selection

%Extraction of the window with the least amount of interframe motion

%Intensity was chosen according to (18) from the paper
intensity=mean(segment,2);

%The extraction follows formula (17) from the paper, adapted due to the
%different sampling frequency in order to have a 25 seconds window 
sum_int=0;
for k=1:length(intensity)-(number_sample+1)
    for i=1:number_sample
        sum_int=sum_int+abs(intensity(k+i)-intensity(k+i+1));  
    end
    interframe(k)=sum_int;
    sum_int=0;
end

[interframe_min,interframe_pos]=min(interframe);

%Extraction of the window for both the RGB signal and the reference signal

segment=segment(interframe_pos:interframe_pos+number_sample,:);
segment_ref=segment_ref(interframe_pos:interframe_pos+number_sample,:);

%extraction of R, G and B values from the final window
R=segment(:,1);
G=segment(:,2);
B=segment(:,3);
%% reference
%The reference of the cardiac frequency was obtained from the finger PPG signal 
%recorded through the ProComp Infiniti encoder, resampled at 115 Hz

freq=(0:1/size(segment,1):1-1/size(segment,1))*fs;
figure('Name','Reference spectrum plot')
reference_spectrum=abs(fft(segment_ref-mean(segment_ref)));
reference_half=reference_spectrum(1:length(reference_spectrum)/2);

%The reference was found in two different ways:
%maximum peak from the spectrum of the signal, in the same way as the
%algorithms
[cardiac_fr_freq_amplitude,spectrum_index]=max(reference_half);
cardiac_fr_freq=freq(spectrum_index);

%computing the peak distance and considering the inverse
[peaks_b,pos_b]=findpeaks(segment_ref,'minpeakdistance',fs*0.5);
cardiac_dist=mean(diff(pos_b))/115;
cardiac_fr_peaks=1/cardiac_dist;

%In order to avoid any bias, we chose to consider, as the actual reference
%throughout the code, the first method
ref_plot=reference_spectrum./max(reference_spectrum);
plot(freq,ref_plot);
title('Reference spectrum');
xlabel('Frequencies (Hz)');
ylabel('Normalized amplitude');
%% RoverG
%First algorithm, application of running average and then processing
%according to (3)

mu=movmean(segment,(fs*2)+1);
normalized=segment./mu;
R_norm=normalized(:,1);
G_norm=normalized(:,2);
B_norm=normalized(:,3);

rover_sgn=G_norm./R_norm-1;

%% XoverY
%Second algorithm (5), we chose to apply a running average as in the previous
%case because it allows for a better reduction of noisy data
X=R-G;
Y=(0.5).*R + (0.5).*G - B;

muX=movmean(X,(fs*2)+1);
muY=movmean(Y,(fs*2)+1);

XoverYsgn=((X./muX)./(Y./muY))-1;

%% fixed

%Third algorithm, since it was not better specified in the paper, we
%applied the running average to the data

R_coeff=0.7682;
G_coeff=0.5121;
B_coeff=0.3841;

Rn=R_norm;
Gn=G_norm;
Bn=B_norm;

%Multiplication of normalized channels to empirical coefficients (7)
Rs=Rn.*R_coeff;
Gs=Gn.*G_coeff;
Bs=Bn.*B_coeff;

%Chrominance signals (9)
Xs=(Rs-Gs)./(R_coeff-G_coeff);
Ys=(Rs+Gs-2.*Bs)/(R_coeff+G_coeff-R_coeff);

%Final formula (8) 
fixed_sgn=(Xs./Ys)-1;

%% XminalfaY
%Fourth algorithm, first the chrominance signals Xs and Ys are bandpassed
filter_opt=fir1(63,[40/(60*(fs/2)) 240/(60*(fs/2))],'bandpass');
Xf=filtfilt(filter_opt,1,Xs);
filter_opt=fir1(63,[40/(60*(fs/2)) 240/(60*(fs/2))],'bandpass');
Yf=filtfilt(filter_opt,1,Ys);

%The disturbances of the final signal were reduced computing alfa
alfa=std(Xf)/std(Yf);
xminalfay_sgn=Xf-alfa.*Yf;


%% ICA
%First benchmark method: the RGB values were initially normalized

R_n_bench=(R-mean(R))./std(R);
G_n_bench=(G-mean(G))./std(G);
B_n_bench=(B-mean(B))./std(B);

%Extracion of 3 independents components
n_comp=3;
rec_ica=rica([R_n_bench G_n_bench B_n_bench],n_comp);
weight_matrix=rec_ica.TransformWeights;
data_matrix=[R_n_bench G_n_bench B_n_bench];

%Reconstruction of components through linear combination between initial
%signal and ICA weights
comp_ica=data_matrix*weight_matrix;
sig_length=length(R_n_bench);
max_sig=0;

%Definition of the component with the highest peak in the normalized
%spectrum, following the indications of the paper
for i=1:3

    spectra_pca=abs(fft(comp_ica(:,i)-mean(comp_ica(:,i))));
    norm_spectra=spectra_pca/sig_length;
    if max(norm_spectra)>max_sig
        max_sig=max(norm_spectra);
        max_ind=i;
    end
end

%Moving average filtering
h=[1/5 1/5 1/5 1/5 1/5];
first_stage=filtfilt(h,1, comp_ica(:,max_ind));



%% PCA
%Second benchmark method

%Initial bandpass filtering
filter_opt_pca_preproc=fir1(127,[40/(60*(fs/2)) 240/(60*(fs/2))],'bandpass',rectwin(128));
R_n_processed=filtfilt(filter_opt_pca_preproc,1,R_n_bench);
G_n_processed=filtfilt(filter_opt_pca_preproc,1,G_n_bench);
B_n_processed=filtfilt(filter_opt_pca_preproc,1,B_n_bench);

weight_pca=pca([R_n_processed G_n_processed B_n_processed]);

%Reconstruction of components through linear combination between initial
%signal and PCA weights
comp_pca=data_matrix*weight_pca;
sig_length=length(R_n_processed);
max_sig_pca=0;

%As before, definition of the component with the highest peak in the normalized
%spectrum

for i=1:3

    spectra_pca=abs(fft(comp_pca(:,i)));
    norm_spectra_pca=spectra_pca/sig_length;
    if max(norm_spectra_pca)>max_sig_pca
        max_sig_pca=max(norm_spectra_pca);
        max_ind_pca=i;
    end

end

%Moving average filtering
h=[1/5 1/5 1/5 1/5 1/5];
first_stage_pca=filtfilt(h,1, comp_pca(:,max_ind_pca));

%% signal plots

%Plot of the signals after processing, these plots are in the time domain

figure('Name','Signals after processing','units','normalized','outerposition',[0 0 1 1])

x_time=(0:1/fs:number_sample/fs);
subplot(3,2,1)
plot(x_time,rover_sgn);
title('RoverG algorithm');
xlabel('Time (s)');
ylabel('Chrominance signal amplitude');
subplot(3,2,2)
plot(x_time,XoverYsgn);
title('XoverY algorithm');
xlabel('Time (s)');
ylabel('Chrominance signal amplitude');
subplot(3,2,3)
plot(x_time,fixed_sgn);
title('Fixed algorithm');
xlabel('Time (s)');
ylabel('Chrominance signal amplitude');
subplot(3,2,4)
plot(x_time,xminalfay_sgn);
title('XminalfaY algorithm');
xlabel('Time (s)');
ylabel('Chrominance signal amplitude');
subplot(3,2,5)
plot(x_time,comp_ica(:,max_ind));
title('ICA algorithm');
xlabel('Time (s)');
ylabel('Signal amplitude');
subplot(3,2,6)
plot(x_time,comp_pca(:,max_ind_pca));
title('PCA algorithm');
xlabel('Time (s)');
ylabel('Signal amplitude');

%% results

%Plots of the spectra. See function 'fun_processing' for more details
fft_coeff=4096;
freq_alg=(0:1/fft_coeff:1-1/fft_coeff)*fs;
freq_alg_halved=freq_alg(1:length(freq_alg)/2);
freq_benchmark_halved=freq(1:length(freq)/2);

figure('Name','Signals spectra','units','normalized','outerposition',[0 0 1 1])

subplot(3,2,1)
[sgn_filtered_rover,fft_coeff_1]=fun_processing(rover_sgn,fs);
sgn_filtered_rover=sgn_filtered_rover./max(sgn_filtered_rover);
plot(freq_alg_halved,sgn_filtered_rover(1:length(sgn_filtered_rover)/2));
title('RoverG algorithm');
xlabel('Frequencies (Hz)');
ylabel('Chrominance signal amplitude');

subplot(3,2,2)
[sgn_filtered_XoverY,fft_coeff_2]=fun_processing(XoverYsgn,fs);
sgn_filtered_XoverY=sgn_filtered_XoverY./max(sgn_filtered_XoverY);
plot(freq_alg_halved,sgn_filtered_XoverY(1:length(sgn_filtered_XoverY)/2));
title('XoverY algorithm');
xlabel('Frequencies (Hz)');
ylabel('Chrominance signal amplitude');

subplot(3,2,3)
[sgn_filtered_fixed,fft_coeff_3]=fun_processing(fixed_sgn,fs);
sgn_filtered_fixed=sgn_filtered_fixed./max(sgn_filtered_fixed);
plot(freq_alg_halved,sgn_filtered_fixed(1:length(sgn_filtered_fixed)/2));
title('Fixed algorithm');
xlabel('Frequencies (Hz)');
ylabel('Chrominance signal amplitude');

subplot(3,2,4)
[sgn_filtered_xminalfay,fft_coeff_4]=fun_processing(xminalfay_sgn,fs);
sgn_filtered_xminalfay=sgn_filtered_xminalfay./max(sgn_filtered_xminalfay);
plot(freq_alg_halved,sgn_filtered_xminalfay(1:length(sgn_filtered_xminalfay)/2));
title('XminalfaY algorithm');
xlabel('Frequencies (Hz)');
ylabel('Chrominance signal amplitude');


subplot(3,2,5)

[ica_spec,fft_coeff_5]=fun_processing(first_stage,fs);
ica_spec=ica_spec./max(ica_spec);
plot(freq_alg_halved,ica_spec(1:length(ica_spec)/2));
title('ICA algorithm');
xlabel('Frequencies (Hz)');
ylabel('Chrominance signal amplitude');

subplot(3,2,6)

[pca_spec,fft_coeff_6]=fun_processing(first_stage_pca,fs);
pca_spec=pca_spec./max(pca_spec);
plot(freq_alg_halved,pca_spec(1:length(pca_spec)/2));
title('PCA algorithm');
xlabel('Frequencies (Hz)');
ylabel('Chrominance signal amplitude');


%% computing snr

%SNR computation according to function 'find_snr', that implements formula (19)
fft_coeff=4096;

snr_rover=find_snr(sgn_filtered_rover,cardiac_fr_freq,fft_coeff,fs);
snr_XoverY=find_snr(sgn_filtered_XoverY,cardiac_fr_freq,fft_coeff,fs);
snr_fixed=find_snr(sgn_filtered_fixed,cardiac_fr_freq,fft_coeff,fs);
snr_xminalfay=find_snr(sgn_filtered_xminalfay,cardiac_fr_freq,fft_coeff,fs);
fft_coeff=length(ica_spec);
snr_ica=find_snr(ica_spec,cardiac_fr_freq,fft_coeff,fs);
fft_coeff=length(pca_spec);
snr_pca=find_snr(pca_spec,cardiac_fr_freq,fft_coeff,fs);

snr_vect=[snr_rover; snr_XoverY; snr_fixed; snr_xminalfay; snr_ica; snr_pca];

%SNR bar plot

figure('Name','SNR')
X = categorical({'RoverG','XoverY','Fixed','XminalfaY','ICA','PCA'});
X = reordercats(X,{'RoverG','XoverY','Fixed','XminalfaY','ICA','PCA'});
bar(X,snr_vect);
ylabel('SNR (dB)');
title('SNR of various methods');
%% find maximum peak from spectra signals

%Computation of cardiac frequency from signals spectra, the highest peak in
%the first half of the spectra is considered as actual cardiac frequency

freq_alg=(0:1/size(sgn_filtered_rover,1):1-1/size(sgn_filtered_rover,1))*fs;

rover_half=sgn_filtered_rover(1:length(sgn_filtered_rover)/2);
[~,rover_index]=max(rover_half);
r_rover=freq_alg(rover_index);

XoverYsgn_half=sgn_filtered_XoverY(1:length(sgn_filtered_XoverY)/2);
[~,XoverYsgn_index]=max(XoverYsgn_half);
r_XoverY=freq_alg(XoverYsgn_index);

fixed_sgn_1_half=sgn_filtered_fixed(1:length(sgn_filtered_fixed)/2);
[~,fixed_sgn_1_index]=max(fixed_sgn_1_half);
r_fixed=freq_alg(fixed_sgn_1_index);

xminalfay_sgn_half=sgn_filtered_xminalfay(1:length(sgn_filtered_xminalfay)/2);
[~,xminalfay_sgn_index]=max(xminalfay_sgn_half);
r_xminalfay=freq_alg(xminalfay_sgn_index);

ica_spec_half=ica_spec(1:length(ica_spec)/2);
[~,ica_spec_index]=max(ica_spec_half);
r_ica=freq_alg(ica_spec_index);

pca_spec_half=pca_spec(1:length(pca_spec)/2);
[ampl,pca_spec_index]=max(pca_spec_half);
r_pca=freq_alg(pca_spec_index);

vect_string=["Reference","RoverG","XoverY","Fixed","XminalfaY","ICA","PCA"];
r_algorithm=[cardiac_fr_freq; r_rover; r_XoverY; r_fixed; r_xminalfay; r_ica; r_pca];

%%
%printing of final values for patient and segment chosen at the beginning
%after conversion from Hz to bpm
result="";
for i = 1:length(vect_string)
    result = sprintf("%s%s: %d BPM\n",result, vect_string(i),round(r_algorithm(i)*60));
end
disp(result);
