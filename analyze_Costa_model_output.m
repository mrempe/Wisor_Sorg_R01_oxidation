function analyze_Costa_model_output(data_file)
% USAGE: analyze_Costa_model_output(data_file)
%
% where data_file is a .mat file containing the output
% of the function run_costa_model_and_save_output.m 
% 

% first load the "data": the output of the simulation
disp('Loading the data file....')
load(data_file);
disp('Done loading!')


% Determine state from the filename
if ~isempty(strfind(data_file,'Wake')) | ~isempty(strfind(data_file,'wake'))
	state = 'Wake';
elseif ~isempty(strfind(data_file,'SWS'))
	state = 'SWS';
elseif ~isempty(strfind(data_file,'REMS')) | ~isempty(strfind(data_file,'REM'))
	state = 'REMS';
else
	error('Unrecognized sleep state in file name');
end



% First filter the data so frequencies below 0.5 Hz are filtered out (High-pass)
load('filter.mat','HP');   					% this is the first HighPass filter FIR I designed using Filter Desing and Analysis App
load('butterworth_HP_filter','HP_butter')	% Butterworth high pass filter
Y             = Y(:,6001:end);  			% First remove the transient from all signals (first 6 seconds)
t 			  = t(6001:end);
data 		  = Y(8,:);          			% Y(8,:) is V_p, a proxy for EEG signal
data_filtered = filter(HP_butter,data);		 % SLOW if using FIR filter.  Try butterworth too. filtering the data using filter.m detrends it.  Just FYI
   
% remove first 6 seconds of data because of edge effect introduced by filtering
data_filtered = data_filtered(20000:end);
t 			  = t(20000:end);
Y 			  = Y(:,20000:end);




% now compute the power in various frequency bands
epoch_length = 10;  % in seconds
freq    	 = 1000/dt;  %10000;  % in Hz
N            = fix(freq*epoch_length);  % number of data points in each epoch
raw_EEG = data_filtered;
raw_EEG_truncated = raw_EEG(1:N*floor(numel(raw_EEG)/N));
raw_EEG_reshaped  = reshape(raw_EEG_truncated,N,length(raw_EEG_truncated)/N);
[spectra_new,freqs_new] = pwelch(raw_EEG_reshaped,hamming(floor(N/2)),[],[],freq);
		
deltaIdx    = (freqs_new>=1  & freqs_new<=4);		%find the indices first
thetaIdx    = (freqs_new>=5  & freqs_new<=9);
gammaIdx    = (freqs_new>=80 & freqs_new<=100);
alphaIdx    = (freqs_new>=8  & freqs_new<=15);
epoch = struct;

epoch.delta = sum(spectra_new(deltaIdx,:),1);   %creates a row vector of delta power in each epoch
epoch.theta = sum(spectra_new(thetaIdx,:),1);
epoch.gamma = sum(spectra_new(gammaIdx,:),1);
epoch.alpha = sum(spectra_new(alphaIdx,:),1);


mean(epoch.delta)
whos Y
% now plot results
figure
plot(t,Y(4,:),'g',t,Y(5,:),'r',t,Y(6,:),'b',t,Y(7,:),'MarkerFaceColor',[1 0.4 0])  % C_E, C_G, C_A, and h 

figure
plot(t,data_filtered)   % V_p, a proxy for the EEG signal.  
title(['Model output EEG signal ', state])
ax=gca;
%ax.YTick = [-90 -60 -30];

figure
plot(epoch.delta)
hold on
plot(epoch.gamma,'g')
plot(epoch.alpha,'r')
legend('delta power','gamma power','alpha power')
hold off
title(state)



% Now make periodogram, overlaying for each sleep state
%data = raw_EEG;
% --- using periodogram       ---
% --- first with all the data ---

window = hamming(length(data_filtered));
[pxxAll,fxAll] = periodogram(data_filtered,window,length(data_filtered),freq,'power');
figure
plot(fxAll,pxxAll)
title(['Simulation data all states using periodogram.  Should be all ', state])
ax=gca;
ax.XLim = [0 100];




% % SWS only 
% figure
% S_indices = find(strcmp(state(:),'S'));
% window = hamming(length(S_indices));
% [pxxS,fxS] = periodogram(data(S_indices),window,length(S_indices),freq,'power');
% plotS=plot(fxS,pxxS,'r');
% ax=gca;
% ax.XLim = [0 100];
% %title('Costa output SWS using periodogram')
% hold on 

% % REMS only 
% R_indices = find(strcmp(state(:),'R'));
% window = hamming(length(R_indices));
% [pxxR,fxR] = periodogram(data(R_indices),window,length(R_indices),freq,'power');
% plotR=plot(fxR,pxxR,'Color',[1 0.5 0]);
% ax=gca;
% ax.XLim = [0 100];
% %title('Costa output REMS using periodogram')
% plotR.Color(4) = 0.2;

% % Wake only 
% W_indices = find(strcmp(state(:),'W'));
% window = hamming(length(W_indices));
% [pxxW,fxW] = periodogram(data(W_indices),window,length(W_indices),freq,'power');
% plotW=plot(fxW,pxxW);
% ax=gca;
% ax.XLim = [0 100];
% title('Experimental data using periodogram')
% hold off
% legend('SWS','REMS','Wake')
% plotW.Color(4) = 0.2;  