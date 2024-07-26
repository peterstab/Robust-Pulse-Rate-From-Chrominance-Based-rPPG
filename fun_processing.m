function [mod_non_norm,fft_coeff]= fun_processing(signal,fs)

%This function is used in order to obtain spectra of the signals after a
%final processing. 
freq=(0:1/size(signal,1):1-1/size(signal,1))*fs;

%Bandpass filter
filter_opt=fir1(127,[40/(60*(fs/2)) 240/(60*(fs/2))],'bandpass',rectwin(128));
signal_filtered=filtfilt(filter_opt,1,signal);

fft_coeff=4096;
%Fft with 4096 points (starting signal was 2875 samples), Hanning window
win=hanning(length(size(signal,1)));
signal_filtered_2=signal_filtered.*win';
yf1=fft(signal_filtered_2-mean(signal_filtered_2),fft_coeff);
mod_non_norm=abs(yf1);


end