function [snr]=find_snr(sgn_filtered, reference,fft_coeff,fs)

%Signal quantification around fundamental frequency (according to reference
%signal) and first harmonic
sgn_filtered=sgn_filtered./max(sgn_filtered);
fundamental=reference;
first_harmonic=fundamental*2;

%bin number was adapted from the paper fft coefficients and sampling frequency ratio (512/20) to the
%signal actual ratio (4096/115). We chose a proportional number
%of bins around the fundamental frequency and the first harmonic comparing
%the two ratios.

fund_bin_index=round((fft_coeff/fs)*fundamental)-3:round((fft_coeff/fs)*fundamental)+3;
first_harm_bin_index=round((fft_coeff/fs)*first_harmonic)-6:round((fft_coeff/fs)*first_harmonic)+7;
signal_index=[fund_bin_index first_harm_bin_index];
k=1;
%binary vector: 1 if the corresponding frequency belongs to the signal
%spectrum, 0 if it is considered as noise
for i=1:fft_coeff
    if i==signal_index(k)
        if k~=length(signal_index)
        k=k+1;
        end
        bin_vect(i)=1;
    else bin_vect(i)=0;
        
    end
end
%binary vector for noise computation, 1 if the bin corresponds to noise, 0
%if it corresponds to signal
bin_vect_snr=1-bin_vect;
size_sgn=size(sgn_filtered,1);
size_bin=size(bin_vect,1);

%check in order to assess correct dimensionality of the vectors
if size_sgn~=size_bin
    sgn_filtered=sgn_filtered';
end

snr_s=sum((sgn_filtered.*bin_vect).^2);

%the spectrum must be considered between 30 and 240 bpm according to (19)
bin_vect_snr=bin_vect_snr(round((fft_coeff/fs)*(30/60)):round((fft_coeff/fs)*(240/60)));
non_norm_temp=sgn_filtered(round((fft_coeff/fs)*(30/60)):round((fft_coeff/fs)*(240/60)));

snr_n=sum((non_norm_temp.*bin_vect_snr).^2);

%final snr
snr=10*log(snr_s/snr_n);