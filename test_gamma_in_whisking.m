% script to read in whisking data and plot gamma vs time and possibly delta vs time
%

% add path for importdatafile.m function
addpath 'C:\Users\wisorlab\Documents\MATLAB\mrempe\Epoch-Based-Processing\Timed_Intervals\export\AutomatedScoring\AutoScore_with_GUI'
addpath '\\FS1\WisorData\MATLAB\Wisor\EpochBasedProcessing\Timed\Excel Export\SleepReport_AnyState_AnyStartTime_1or2EEGs\'  % for rescoreQuietWakeActiveWake 
%
directory = '\\FS1\WisorData\Jonathan Data\LactateWhisker-FftsAndEdfs\';
files = dir(strcat(directory,'*.txt'));


EmgColumn = 42;

for i = 1: length(files)
	[data,textdata] = importdatafile([directory files(i).name]);

	Gamma1{i} 	  = data(:,44);
	Gamma2{i} 	  = data(:,46);
	SleepState{i} = RescoreQuietVsActiveWake(char(textdata(:,2)),data(:,EmgColumn),0.33,0.66,i,files);
	leng(i) = length(Gamma1{i});

end


% set up start and end of tickle sessions so we can plot those:
tickle_start = [8649 9009 9369 9729 10089 10449];
tickle_end = tickle_start + 180;


% average gamma across recordings (may need to cut all recordings to be length of shortest recording)
% and compute SEM  
% find shortest recording and truncate all recordings to this length
min_length = min(leng);
for i=1:length(files)
	Gamma1{i}     = Gamma1{i}(1:min_length);
	Gamma2{i}     = Gamma2{i}(1:min_length);
	SleepState{i} = SleepState{i}(1:min_length);
end

% make a big matrix for Gamma1, Gamma2
Gamma1_matrix = Gamma1{1};
Gamma2_matrix = Gamma2{1};
for i=2:length(Gamma1)
	Gamma1_matrix = [Gamma1_matrix Gamma1{i}];
	Gamma2_matrix = [Gamma2_matrix Gamma2{i}];
end 


average_gamma1 = mean(Gamma1_matrix,2);
average_gamma2 = mean(Gamma2_matrix,2);

dt=10/60/60;  %epoch length in hours
timevec = 0:dt:dt*length(average_gamma1);
timevec = timevec(1:end-1);
% plot gamma vs time Also color code by sleep state
figure
plot(average_gamma1)

% add rectangles to show when tickling happened
hold on
for i=1:length(tickle_start)
	p = patch([tickle_start(i) tickle_end(i) tickle_end(i) tickle_start(i)],[0 0 350 350],[0.5725 0.2863 0]);
	p.EdgeColor='none';
	p.FaceAlpha = 0.7;
end 

hold off
legend('Gamma power averaged across recordings')
title ('Gamma power (80-90 Hz) in channel 1')

figure
plot(average_gamma2)
hold on
for i=1:length(tickle_start)
	p = patch([tickle_start(i) tickle_end(i) tickle_end(i) tickle_start(i)],[0 0 200 200],[0.5725 0.2863 0]);
	p.EdgeColor='none';
	p.FaceAlpha = 0.7;
end 
hold off
legend('Gamma power averaged across recordings')
title ('Gamma power (80-90 Hz) in channel 2')



% plot one gamma trace with color for sleep state
num=11;
figure
hold on 
for i=1:size(Gamma1_matrix,1)-1
	if strcmp(SleepState{num}(i),'W')
		color = 'r';
	elseif strcmp(SleepState{num}(i),'A')
		color = [0.85 0.43 0];
	elseif strcmp(SleepState{num}(i),'Q')
		color = [0.4 0.7 1];
	elseif strcmp(SleepState{num}(i),'S')
		color = 'g';
	elseif strcmp(SleepState{num}(i),'R')
		color = 'b';
	else
		color = 'k';
	end 

	L=line([timevec(i) timevec(i+1)],[Gamma1_matrix(i,num) Gamma1_matrix(i+1,num)]);
	L.Color = color;
end
for i=1:length(tickle_start)
	p = patch([timevec(tickle_start(i)) timevec(tickle_end(i)) timevec(tickle_end(i)) timevec(tickle_start(i))],[0 0 200 200],[0.5725 0.2863 0]);
	p.EdgeColor='none';
	p.FaceAlpha = 0.7;
end 
hold off


% compute averages during recovery sessions?   