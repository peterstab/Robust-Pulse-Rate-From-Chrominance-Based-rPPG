function [final_results]= computing_data(data_rgb,reference)

%computation of all cardiac frequencies and indexes for every patient,
%phase and algorithm
wait = size(data_rgb,2)-1;
h = waitbar(0,'Please wait...');
for i=1:size(data_rgb,2)-1
    
    
    
    %extraction of each phase from the starting RGB signals
    
    data_tot=data_rgb{i}(end-65549:end,:);
    bvp_temp=reference(2:end,i);
    fs=115;
    indexes=[fs*60*5 fs*60*7 fs*60*8 fs*60*8.5 fs*60*9 fs*60*9.5];
    start_baseline=data_tot(1:indexes(1),:);
    paced_breathing=data_tot(indexes(1)+1:indexes(2),:);
    rest1=data_tot(indexes(2)+1:indexes(3),:);
    hand_grip1=data_tot(indexes(3)+1:indexes(4),:);
    rest2=data_tot(indexes(4)+1:indexes(5),:);
    hand_grip2=data_tot(indexes(5)+1:indexes(6),:);
    
    %extraction of each phase from the starting reference signals
    start_baseline_ref=bvp_temp(1:indexes(1),:);
    paced_breathing_ref=bvp_temp(indexes(1)+1:indexes(2),:);
    rest1_ref=bvp_temp(indexes(2)+1:indexes(3),:);
    hand_grip1_ref=bvp_temp(indexes(3)+1:indexes(4),:);
    rest2_ref=bvp_temp(indexes(4)+1:indexes(5),:);
    hand_grip2_ref=bvp_temp(indexes(5)+1:indexes(6),:);
    number_sample=2875;
    j=1;
    
    %application of 'segment_analysis' function to extract every cardiac
    %frequency and every snr for every phase
    
    %result_matrix will be a tensor with 
    %first dimension: phase
    %second dimension: algorithm (the first column is reserved to the
    %reference values)
    %third dimension: patient
    
    [freq_ref,algorithm_vect,snr_seg]=segment_analysis(start_baseline,start_baseline_ref,number_sample);
    result_matrix(j,:,i)=[freq_ref; algorithm_vect];
    snr_matrix(j,:,i)=snr_seg;
    j=j+1;
    [freq_ref,algorithm_vect,snr_seg]=segment_analysis(paced_breathing,paced_breathing_ref,number_sample);
    result_matrix(j,:,i)=[freq_ref; algorithm_vect];
    snr_matrix(j,:,i)=snr_seg;
    j=j+1;
    [freq_ref,algorithm_vect,snr_seg]=segment_analysis(rest1,rest1_ref,number_sample);
    result_matrix(j,:,i)=[freq_ref; algorithm_vect];
    snr_matrix(j,:,i)=snr_seg;
    j=j+1;
    [freq_ref,algorithm_vect,snr_seg]=segment_analysis(hand_grip1,hand_grip1_ref,number_sample);
    result_matrix(j,:,i)=[freq_ref; algorithm_vect];
    snr_matrix(j,:,i)=snr_seg;
    j=j+1;
    [freq_ref,algorithm_vect,snr_seg]=segment_analysis(rest2,rest2_ref,number_sample);
    result_matrix(j,:,i)=[freq_ref; algorithm_vect];
    snr_matrix(j,:,i)=snr_seg;
    j=j+1;
    [freq_ref,algorithm_vect,snr_seg]=segment_analysis(hand_grip2,hand_grip2_ref,number_sample);
    result_matrix(j,:,i)=[freq_ref; algorithm_vect];
    snr_matrix(j,:,i)=snr_seg;
    waitbar(i/wait,h);
   
   

end

%% RMSE and STD

%RMSE and STD computation from the result matrix obtained above
%For each index, the result is a 6x6 matrix, each row representing a phase
%and each column an algorithm

rmse_value=zeros(6,6);
for pat=1:size(result_matrix,3)
    for method=2:size(result_matrix,2)
        for segment=1:size(result_matrix,1)
rmse_value_temp=(result_matrix(segment,1,pat)*60-result_matrix(segment,method,pat)*60)^2;
rmse_value(segment,method-1)=rmse_value(segment,method-1)+rmse_value_temp;
std_tensor(segment, method-1, pat)=abs(result_matrix(segment,1,pat)*60-result_matrix(segment,method,pat)*60);
% std_tensor(segment, method-1, pat)=result_matrix(segment,method,pat);
        end
    end
end
rmse_final=sqrt((rmse_value./size(result_matrix,3)));
std_final=std(std_tensor, [], 3);

%% Pearson and Slope coefficients

%Used functions: corrcoef with, as input, the cardiac frequencies obtained
%from a specific algorithm and a specific phase for every patient (total
%length: 20) compared to the cardiac frequences of the reference for every
%patient (same length of the previous vector)
%Pearson matrix will be a 6x6 matrix, each row will represent a phase, each
%column an algorithm

%The same procedure was followed for the computation of the slopes, using
%polyfit function and finding the slope coefficient

pearson_matrix=zeros(6,6);
slope_matrix=zeros(6,6);
for method=2:size(result_matrix,2)
    for segment=1:size(result_matrix,1)
        pearson_temp=corrcoef(squeeze(result_matrix(segment,1,:)),squeeze(result_matrix(segment,method,:)));
        pearson_matrix(segment, method-1)=pearson_temp(1,2);
       slope_temp=polyfit(squeeze(result_matrix(segment,1,:)),squeeze(result_matrix(segment,method,:)),1);
       slope_matrix(segment, method-1)=slope_temp(1);
    end
end

close(h)

%Final struct containing cardiac frequencies and indexes

final_results=struct('cardiac_frequencies',result_matrix,'SNR',snr_matrix,'RMSE',rmse_final,'STD',std_final,'Pearson_coeff',pearson_matrix,'slope_of_linear_fit',slope_matrix);
end