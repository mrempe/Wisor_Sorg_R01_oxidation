function [S,state,long_wake_episode_timings,sleep_measure_averages,sleep_measure_stds] = two_process_model_with_markov_chain(total_time,input_params,shift,makeplots)
% USAGE: [S,state,long_wake_episode_timings,sleep_measure_averages,sleep_measure_stds] = two_process_model_with_markov_chain(total_time,input_params,shift,makeplots)
%
% This function simulates the changes in sleep state using a markov chain model similar to Kemp and Kamphuisen SLEEP 1986
% The markov chain generates sleep states based on the values of a homeostat and the Circadian curve.  
% The homeostat goes up with time constant taui when the sleep state is W or R and goes down with time constant taud when the sleep state is S.
% 
% This was originally developed to model rat sleep in a simulated shift work protocol (published in Rempe et al Neurobio Sleep Circ Rhy 2018)
% 
% for the AW (active worker) case:
%  [S,num_wake_episodes,mean_SWS_length,mean_REM_length]=two_process_model_with_markov_chain(134,8.6,3.2,'AW',1)	
%
% for the RW (resting worker) case:
%  [S,num_wake_episodes,mean_SWS_length,mean_REM_length]=two_process_model_with_markov_chain(134,8.6,3.2,'RW',1)
%
%


if nargin == 3
	sleep_dep_start_stop_times = [];
	makeplots = 0;
end 



taui_baseline = input_params.taui_baseline;
taud_baseline = input_params.taud_baseline;
taui_work     = input_params.taui_work;
taud_work     = input_params.taud_work;
 
alertness_duration_scale_factor_B       = input_params.alertness_duration_scale_factor_B;
sleepiness_duration_scale_factor_SWS_B  = input_params.sleepiness_duration_scale_factor_SWS_B;  
sleepiness_duration_scale_factor_REMS_B = input_params.sleepiness_duration_scale_factor_REMS_B; 
alertness_duration_scale_factor_W       = input_params.alertness_duration_scale_factor_W;
sleepiness_duration_scale_factor_SWS_W  = input_params.sleepiness_duration_scale_factor_SWS_W;  
sleepiness_duration_scale_factor_REMS_W = input_params.sleepiness_duration_scale_factor_REMS_W;  




% --- Set up the time vector (with increments of 10 seconds)
t=0:10/60/60:total_time;
dt=10/60/60;
phase = mod(t,24);
num_days = total_time/24;
% -----------------------------------------------------------


% --- Set up the sleep deprivation, if applied
if strcmp(shift,'AW')
	sleep_dep_start_stop_times = [38 46 62 70 86 94 110 118];
elseif strcmp(shift,'RW')
	sleep_dep_start_stop_times = [26 34 50 58 74 82 98 106];
else
	sleep_dep_start_stop_times = [];
end 
sleep_dep = zeros(size(t));
if ~isempty(sleep_dep_start_stop_times)
	if mod(length(sleep_dep_start_stop_times),2) ~= 0
		error('sleep_dep_start_stop_times must have an even number of elements')
	end 

	sleep_dep_start_times = sleep_dep_start_stop_times(1:2:end-1);
	sleep_dep_end_times	  = sleep_dep_start_stop_times(2:2:end);

	for i=1:length(sleep_dep_start_times)
		ind_sleep_dep_start(i) = find(abs(t-sleep_dep_start_times(i))<1e-6);
		ind_sleep_dep_end(i)   = find(abs(t-sleep_dep_end_times(i))<1e-6);
		sleep_dep(ind_sleep_dep_start(i):ind_sleep_dep_end(i)) = 1;  % 23-25 hours
	end 

end
% -----------------------------------------------------------



% -----  Circadian curve -----
circ_amp = 0.25; %0.128;
circ_curve = circ_amp*sin((pi/12)*-t);
C = circ_curve;
% ----------------------------



num_simulations = 50; %50;  % was 100
% ----- Initialize S, state and C
state(1,:) = char(repmat('S',1,num_simulations));   % initialize the first row of state matrix to S
S=zeros(length(t),num_simulations);
S(1,:) = 0.3;					% starting value for S, the homeostat
% -------------------------------


% Set up the random variables in a power law for W episode duration
% I found these values using the MATLAB function dfittool.m. Units are in seconds.
pd_w = makedist('GeneralizedPareto','k',0.231707,'sigma',217.976,'theta',1);
% or try k=0.231707 and sigma = 217.976 k=0.267395, sigma=20.1499
pd_s = makedist('Burr','alpha',177.113,'c',4.10314,'k',1.48324);
pd_r = makedist('Burr','alpha',172.817,'c',2.629,'k',4.75);



taui = taui_baseline*ones(size(t));
taud = taud_baseline*ones(size(t));
alertness_duration_scale_factor       = ones(size(t));
sleepiness_duration_scale_factor_SWS  = ones(size(t));
sleepiness_duration_scale_factor_REMS = ones(size(t));


taui(1:8640)   = taui_baseline;
taui(8641:end) = taui_work;
taud(1:8640)   = taud_baseline;
taud(8641:end) = taud_work;



alertness_duration_scale_factor(1:8640) 		= alertness_duration_scale_factor_B; 
alertness_duration_scale_factor(8641:end)       = alertness_duration_scale_factor_W;
sleepiness_duration_scale_factor_SWS(1:8640)    = sleepiness_duration_scale_factor_SWS_B;
sleepiness_duration_scale_factor_SWS(8641:end)  = sleepiness_duration_scale_factor_SWS_W;
sleepiness_duration_scale_factor_REMS(1:8640)   = sleepiness_duration_scale_factor_REMS_B;
sleepiness_duration_scale_factor_REMS(8641:end) = sleepiness_duration_scale_factor_REMS_W; 





long_wake_episode_timings = zeros(size(t));
% Run the MC simulation for several different runs, each column corresponds to a different run
for run = 1:num_simulations
	num_R_runs = 0;
	i=1;
	long_wake_counter(run) = 0;
	long_wake_episode_timings = zeros(size(t));
	while i<length(t)
		
		
		sleepiness = S(i,run)-C(i);
		alertness = 1-S(i,run)+C(i);
		

		if state(i,run)=='W'
			% rare long wake episode (uniformly distributed between 5 minutes and 30 minutes which is 30 epochs and 180 epochs)  (b-a)*u+a gives
			% a uniform random variable on [a b].  
			
			if alertness > (8-0.5)*rand + 0.5
				W_SW = round((180-30)*rand + 30);  % a uniform random variable between 30 and 180 epochs 
				Wshortest_epochs = W_SW;
				if i+Wshortest_epochs > length(t) 
					Wshortest_epochs = length(t)-i; 
				end % truncate if the wait time puts you past the end of t
				next_state = 'S';
				state(i+1:i+W_SW,run)='W';  
				state(i+W_SW+1,run)='S';     % after long run of Wake is finished set the next state to S
				for j=1:W_SW				 % update the homeostat
					S(i+j,run)=1-(1-S(i+j-1,run))*exp(-dt/taui(i));
				end
				long_wake_counter(run) = long_wake_counter(run)+1;
				long_wake_episode_timings(i) = 1;    
				
			else  % normal markov chain 
				
				% choose the next state and the bout duration separately
				if rand < 0.135*sleepiness-0.0135
					% go to REMS
					next_state = 'R';
				else
					next_state = 'S';
				end 
				Wshortest_epochs = round((random(pd_w)/10)*alertness_duration_scale_factor(i)*(alertness+0.3)-1);  % scaling by alertness since W bout durations are a function of alertness
				if Wshortest_epochs <= 0 
					Wshortest_epochs = 1;
					%warning('Wshortest_epochs<=0 in W') 
				end 
				if i+Wshortest_epochs > length(t) 
					Wshortest_epochs = length(t)-i; 
				
				end % truncate if the wait time puts you past the end of t
				
				state(i+1:i+Wshortest_epochs,run)='W';
				state(i+Wshortest_epochs+1,run)=next_state;  % after the run of Wake is finished set the next state depending on which W was shortest
				for j=1:Wshortest_epochs				     % update the homeostat
					S(i+j,run) = 1-(1-S(i+j-1,run))*exp(-dt/taui(i));
				end
				
			end % end of normal markov chain case (not long W episode)


		elseif state(i,run)=='S'
			% choose the next state and the bout duration separately
			% looking at the baseline data, S is initally equally likely to transition to R or W, but throughout
			% the 24 hour baseline, W becomes more likely in a linear fashion.  
			% the line is given by y = (-0.17/24)*(t) + 0.5
			if rand < -0.007*phase+0.51
				next_state = 'R';
			else
				next_state = 'W';
			end 
			Wshortest_epochs = round((random(pd_s)/10)*sleepiness_duration_scale_factor_SWS(i)*(sleepiness+0.7)-1); % the /10 is to get epochs not seconds.  multiplty by sleepiness since S bout duration scales with sleepiness. 0.7 to get sleepiness to average 1.
																	% may need a scale factor here: scale*(sleepiness+0.7)
																	% minus one because we already have one epoch of S from end of previous run
			if Wshortest_epochs <= 0 
				Wshortest_epochs = 1;
				end 
			if i+Wshortest_epochs > length(t) 
				Wshortest_epochs = length(t)-i; 
			end % truncate if the wait time puts you past the end of t
			
			state(i+1:i+Wshortest_epochs,run)='S';
			state(i+Wshortest_epochs+1,run)=next_state;  % after the run of S is finished set the next state depending on which W was shortest
			for j=1:Wshortest_epochs
				S(i+j,run) = S(i+j-1,run)*exp(-dt/taud(i));
			end
			
	

		elseif state(i,run)=='R'
			num_R_runs = num_R_runs+1;
			
			if rand < 0.5  % R is equally likely to go to W or S (found from baseline data)
				next_state = 'W';
			else
				next_state = 'S';
			end 
			Wshortest_epochs = round((random(pd_r)/10)*sleepiness_duration_scale_factor_REMS(i)*(sleepiness+0.7)-1); %minus 1 because we already have one epoch that is R
			
			if Wshortest_epochs <= 0 
				Wshortest_epochs = 1;
				end 
			if i+Wshortest_epochs > length(t) 
				Wshortest_epochs = length(t)-i; 
			end % truncate if the wait time puts you past the end of t
			
			state(i+1:i+Wshortest_epochs,run) = 'R';
			state(i+Wshortest_epochs+1,run)   = next_state;  % after the run of REMS is finished set the next state depending on which W was shortest
			for j=1:Wshortest_epochs				 		 % update the homeostat
				S(i+j,run) = 1-(1-S(i+j-1,run))*exp(-dt/taui(i));
			end

		end
		if strcmp(next_state,'W') || strcmp(next_state,'R')
			S(i+Wshortest_epochs+1,run) = 1-(1-S(i+Wshortest_epochs,run))*exp(-dt/taui(i));
		else
			S(i+Wshortest_epochs+1,run) = S(i+Wshortest_epochs,run)*exp(-dt/taud(i)); 
		end

		% Finally, override the Markov Chain sleep state if sleep deprivation is turned on
		if i+Wshortest_epochs+1 < length(t)
			if sum(find(sleep_dep(i+1:i+Wshortest_epochs+1)))>0
				sleep_dep_locs_this_cycle = find(sleep_dep(i+1:i+Wshortest_epochs+1));
				for j=sleep_dep_locs_this_cycle
					state(i+1+j-1,run) = 'W';
					S(i+1+j-1,run) = 1-(1-S(i+1+j-2,run))*exp(-dt/taui(i));
				end
			end
		end

		% update sleepiness and alertness values, since S has been updated
		for j=1:Wshortest_epochs				 % update the homeostat
			% new
			sleepiness = S(i+j,run)-C(i+j);
			alertness = 1-S(i+j,run)+C(i);
		end



		i=i+Wshortest_epochs+1;
	
		Wshortest_epochs = [];  % faster than clear Wshortest_epochs, same result
	end

end  % looping over num_simulations


 


for j=1:num_simulations
	wake_percentage(j) = length(find(state(:,j)=='W'))/size(state,1);
	SWS_percentage(j)  = length(find(state(:,j)=='S'))/size(state,1);
	REM_percentage(j)  = length(find(state(:,j)=='R'))/size(state,1);
end

av_wake_percentage = mean(wake_percentage);
av_SWS_percentage  = mean(SWS_percentage);
av_REM_percentage  = mean(REM_percentage);



% ---- Compute average episode length -----
[sleep_measure_averages,sleep_measure_stds] = group_and_plot_Wepisodes_SWS_REMS_duration(shift,state,num_simulations,0);




if makeplots
	
	% ---- first make a shaded plot showing average trace of the homeostat ----
	if length(S) > length(t)
		S = S(1:length(t),:);
		tnew = t;
		circ_curvenew = circ_curve;
	else
		tnew = t;
		circ_curvenew = circ_curve;
	end 

	hfig1=figure;
	err = std(S,0,2);
	set(hfig1,'Position',[200 700 1400 275])
	fill([tnew';flipud(tnew')],[mean(S,2)-err;flipud(mean(S,2)+err)],[0.9 0.9 0.9],'linestyle','none');
	l=line(tnew',mean(S,2));
	set(l,'linewidth',1.5)
	ax=gca;
	ax.XTick = 0:12:t(end);
	set(gca,'XTick',0:24:t(end));	% set major tick marks at 24 hour intervals
	set(gca,'YTick',0:0.1:1);
	ax.XAxis.MinorTick = 'on';
	ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):8:ax.XAxis.Limits(2);  % set minor tick marks at 8 hour intervals
	ax.TickDir = 'out';
	ax.TickLength = [0.02 0.02];
	ax.LineWidth = 1.5;
	ax.FontSize = 16;
	set(gca,'box','off')
	set(gca,'color','none')
	ylabel('Homeostat')

	% add gray bars indicating work periods
	hold on
	if ~isempty(sleep_dep_start_stop_times)
		gray = [0.31 0.31 0.31];
		for i=1:length(sleep_dep_start_times)
			p = patch([sleep_dep_start_times(i) sleep_dep_end_times(i) sleep_dep_end_times(i) sleep_dep_start_times(i)],[0 0 0.8 0.8],gray);
			p.EdgeColor='none';
			p.FaceAlpha = 0.6;
		end 
	end
	rectangle('Position',[0 0.8 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[12 0.8 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[24 0.8 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[36 0.8 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[48 0.8 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[60 0.8 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[72 0.8 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[84 0.8 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[96 0.8 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[108 0.8 12 0.1],'EdgeColor','none','FaceColor','k')
	if t(end)>120
		rectangle('Position',[120 0.8 12 0.1],'EdgeColor','none','FaceColor','y')
	end
	hold off 
	
% ---  Sleepiness and Alertness   ---------------------------------------------------------
	hfig2=figure;
	set(hfig2,'Position',[200 350 1400 275])
	p=plot(tnew',mean(S,2)-circ_curvenew',tnew',1-mean(S,2)+circ_curvenew','r');
	hold on
	longw_locs = find(long_wake_episode_timings);
	plot(t(longw_locs),0.9*long_wake_episode_timings(longw_locs),'.','MarkerSize',16)
	hold off
	set(p,'linewidth',1.5)
	ax=gca;
	ax.XTick = 0:12:t(end);
	h_legend = legend('Sleepiness','Alertness');
	legend(ax,'boxoff')
	set(h_legend,'FontSize',16)
	ax=gca;
	ax.XTick = 0:12:t(end);
	set(gca,'XTick',0:24:t(end));	% set major tick marks at 24 hour intervals
	set(gca,'YTick',0:0.2:1);
	ax.XAxis.MinorTick = 'on';
	ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):8:ax.XAxis.Limits(2);  % set minor tick marks at 8 hour intervals
	ax.TickDir = 'out';
	ax.TickLength = [0.02 0.02];
	ax.LineWidth = 1.5;
	ax.FontSize = 16;
	set(gca,'box','off')
	set(gca,'color','none')
	xlabel('TIME (H)')
	axis([0 140 -0.1 1.75])

	% add gray bars indicating work periods
	hold on
	if ~isempty(sleep_dep_start_stop_times)
		gray = [0.31 0.31 0.31];
		for i=1:length(sleep_dep_start_times)
			p = patch([sleep_dep_start_times(i) sleep_dep_end_times(i) sleep_dep_end_times(i) sleep_dep_start_times(i)],[0 0 1 1],gray);
			p.EdgeColor='none';
			p.FaceAlpha = 0.6;
		end 
	end
	rectangle('Position',[0 1 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[12 1 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[24 1 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[36 1 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[48 1 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[60 1 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[72 1 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[84 1 12 0.1],'EdgeColor','none','FaceColor','k')
	rectangle('Position',[96 1 12 0.1],'EdgeColor','none','FaceColor','y')
	rectangle('Position',[108 1 12 0.1],'EdgeColor','none','FaceColor','k')
	if t(end)>120
		rectangle('Position',[120 1 12 0.1],'EdgeColor','none','FaceColor','y')
	end
	hold off 
	


	% ---- Then make a colored shaded plot of percentages of states, like Figure 5 in the experimental manuscript
	make_shaded_state_percentages_plot(t,state,sleep_dep_start_stop_times)

end

clear long_wake_counter longw_locs

disp('end of two_process_model_with_markov_chain')

 