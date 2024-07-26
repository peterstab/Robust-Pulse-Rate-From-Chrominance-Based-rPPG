function []=computing_results(phase,result_struct)
%Computation of indexes depending on the phase chosen by the user

titleVect = ["Baseline" "Paced breathing" "Rest1" "Handgrip1" "Rest2" "Handgrip2"];


SNR_algorithms=mean(squeeze(result_struct.SNR(phase,:,:)),2);
fig=figure('Name',titleVect(phase));
X = categorical({'RoverG','XoverY','Fixed','XminalfaY','ICA','PCA'});
X = reordercats(X,{'RoverG','XoverY','Fixed','XminalfaY','ICA','PCA'});
coordinates=[0 350; 0 0; 500 350; 500 0; 1000 350; 1000 0];
ax = uiaxes(fig);
ax.Position(2) = 115;
fig.Position(3)=480;
bar(ax,X,SNR_algorithms);
ylabel('SNR (dB)');
title(titleVect(:,phase)+' results');
movegui(fig,coordinates(phase,:));
%%
matrix_temp(1,:)=result_struct.RMSE(phase,:);
matrix_temp(2,:)=result_struct.STD(phase,:);
matrix_temp(3,:)=result_struct.slope_of_linear_fit(phase,:);
matrix_temp(4,:)=result_struct.Pearson_coeff(phase,:);

uit = uitable(fig,'Data',matrix_temp,'Position',[0 0 480 98]);
uit.ColumnName  = {'RoverG','XoverY','Fixed','XminalfaY','ICA','PCA'};
uit.RowName={'RMSE (bpm)','STD (bpm)','SLOPE','PEARSON'};
end
