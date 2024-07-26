%%
%open workspace and compute every cardiac frquence, RMSE, STD, and other
%indexes for every patient, algorithm and phase of the signal
clear all
close all
load('starting_workspace.mat');

final_results=computing_data(total_data,bvp);

%%

%phases and relative indexes


%baseline (5min)          -> phase=1
%paced breathing (2min)   -> phase=2
%rest1 (1min)             -> phase=3
%handgrip1 (30sec)        -> phase=4
%rest2  (30sec)           -> phase=5
%handgrip2 (30sec)        -> phase=6

% or cycle for results of every phase



for phase=1:6
computing_results(phase,final_results);
end


%% Only on selected patients

% We can select which patients can be included in the results. We can choose the 10 patients
% whose RGB are centered around 3 different values and discard other RGB signal
% with similar values over channels.
% Results improve

fs = 115;
time = 0:1/fs:(length(total_data{1,1}(:,1))-1)/fs;

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:20
subplot(4,5,i); 
plot(time, total_data{1,i}(:,1), 'red')
hold on 
plot(time, total_data{1,i}(:,2), 'green')
plot(time, total_data{1,i}(:,3), 'blue')
xlabel('time [s]'); ylabel ('pixel intensity [-]');

title( sprintf('RGB Patient %g', i) )

end
%%
selected_pat = [2 6 7 9 12 14 15 16 17 18 21];

new_data = total_data(1, selected_pat);
new_bvp = bvp(:, selected_pat);
new_final_results=computing_data(new_data,new_bvp);

for phase=1:6
computing_results(phase,new_final_results);
end



