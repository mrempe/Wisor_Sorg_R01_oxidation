% script to load and plot Ashley's astrocyte Ca data 
% format: 
% Timestamp		ttl_start	zeitgeiber	treatment	epochNo	SleepStage	Ca data for ROIs

filename = 'GSTIM-M1_ROI_Ca2_sleepstages.xlsx';

% read in timestamps
[~,timestamps,~] = xlsread(filename,1,'A2:A652');

% read in the sleep state
[~,sleepstate,~] = xlsread(filename,1,'E2:E652');


% read in the data traces
traces = xlsread(filename,1,'G2:FF652');



% First plot data from every stream using one color for each stream
figure
for i = 1:size(traces,2)
	plot(traces(:,i))
	hold on
end




% Next plot (a subset of) the data using colors for sleep state