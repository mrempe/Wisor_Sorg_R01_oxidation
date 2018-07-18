function [averages,stndev] = group_and_plot_Wepisodes_SWS_REMS_duration(shift,state,num_simulations,makeplots)
% USAGE: [averages,stndev] = group_and_plot_Wepisodes_SWS_REMS_duration(shift,state,num_simulations,makeplots)


if strcmp(shift,'AW') || strcmp(shift,'none')
	W1_indices = 16561:22321;  % 16 hours
	W2_indices = 25201:30961;
	W3_indices = 33841:39601;
	W4_indices = 42481:48241;
elseif strcmp(shift,'RW')
	W1_indices = 12241:18001;
	W2_indices = 20881:26641;
	W3_indices = 29521:35281;
	W4_indices = 38161:43921;
else
	error('Shift must either be ''AW'',''RW'', or ''none'' ')
end 




for j=1:num_simulations
	baseline_runs = contiguous(state(1:8641,j),['W' 'S' 'R']);
	W1_runs 	  = contiguous(state(W1_indices,j),['W' 'S' 'R']);
	W2_runs 	  = contiguous(state(W2_indices,j),['W' 'S' 'R']);
	W3_runs 	  = contiguous(state(W3_indices,j),['W' 'S' 'R']);
	W4_runs 	  = contiguous(state(W4_indices,j),['W' 'S' 'R']);

	% only consider W and SWS runs of 3 epochs or longer and REMS runs of 2 epochs or longer
	% baseline_runs = remove_short_runs(baseline_runs);
	% W1_runs = remove_short_runs(W1_runs);
	% W2_runs = remove_short_runs(W2_runs);
	% W3_runs = remove_short_runs(W3_runs);
	% W4_runs = remove_short_runs(W4_runs);



	try 
		B_wake_runs{j} = baseline_runs{1,2};
	catch
		B_wake_runs{j} = [];
	end
	try 
		B_sws_runs{j}  = baseline_runs{2,2};
	catch 
		B_sws_runs{j} = [];
	end
	try
		B_rems_runs{j}  = baseline_runs{3,2};
	catch 
		B_rems_runs{j} = [];
	end
	try 
		W1_wake_runs{j} = W1_runs{1,2};
	catch 
		W1_wake_runs{j} = [];
	end
	try
		W1_sws_runs{j} = W1_runs{2,2};
	catch 
		W1_sws_runs{j} = [];
	end 
	try
		W1_rems_runs{j} = W1_runs{3,2};
	catch
		W1_rems_runs{j} = [];
	end
	try 
		W2_wake_runs{j} = W2_runs{1,2};
	catch
		W2_wake_runs{j} = [];
	end
	
	try
		W2_sws_runs{j} = W2_runs{2,2};
	catch
		W2_sws_runs{j} = [];
	end
	try
		W2_rems_runs{j} = W2_runs{3,2};
	catch
		W2_rems_runs{j} = [];
	end 
	try 
		W3_wake_runs{j} = W3_runs{1,2};
	catch
		W3_wake_runs{j} = [];
	end
	try
		W3_sws_runs{j} = W3_runs{2,2};
	catch
		W3_sws_runs{j} = [];
	end
	try
		W3_rems_runs{j} = W3_runs{3,2};
	catch 
		W3_rems_runs{j} = [];
	end
	try
		W4_wake_runs{j} = W4_runs{1,2};
	catch 
		W4_wake_runs{j} = [];
	end
	try
		W4_sws_runs{j} = W4_runs{2,2};
	catch
		W4_sws_runs{j} = [];
	end
	try 
		W4_rems_runs{j} = W4_runs{3,2};
	catch
		W4_rems_runs{j} = [];
	end

	% Number of wake epochs (not episodes)
	B_num_wake_epochs(j)  = length(find(state(1:8641,j)=='W'));
	W1_num_wake_epochs(j) = length(find(state(W1_indices,j)=='W'));
	W2_num_wake_epochs(j) = length(find(state(W2_indices,j)=='W'));
	W3_num_wake_epochs(j) = length(find(state(W3_indices,j)=='W'));
	W4_num_wake_epochs(j) = length(find(state(W4_indices,j)=='W'));

	% Number of sws epochs 
	B_num_sws_epochs(j)  = length(find(state(1:8641,j)=='S'));
	W1_num_sws_epochs(j) = length(find(state(W1_indices,j)=='S'));
	W2_num_sws_epochs(j) = length(find(state(W2_indices,j)=='S'));
	W3_num_sws_epochs(j) = length(find(state(W3_indices,j)=='S'));
	W4_num_sws_epochs(j) = length(find(state(W4_indices,j)=='S'));
	
	% Number of rems epochs
	B_num_rems_epochs(j)  = length(find(state(1:8641,j)=='R'));
	W1_num_rems_epochs(j) = length(find(state(W1_indices,j)=='R'));
	W2_num_rems_epochs(j) = length(find(state(W2_indices,j)=='R'));
	W3_num_rems_epochs(j) = length(find(state(W3_indices,j)=='R'));
	W4_num_rems_epochs(j) = length(find(state(W4_indices,j)=='R'));


	% Number of wake episodes (not epochs)
	B_num_wake_episodes(j)  = size(B_wake_runs{j},1);
	W1_num_wake_episodes(j) = size(W1_wake_runs{j},1);
	W2_num_wake_episodes(j) = size(W2_wake_runs{j},1);
	W3_num_wake_episodes(j) = size(W3_wake_runs{j},1);
	W4_num_wake_episodes(j) = size(W4_wake_runs{j},1);

	% Number of SWS episodes
	B_num_sws_episodes(j)  = size(B_sws_runs{j},1);
	W1_num_sws_episodes(j) = size(W1_sws_runs{j},1);
	W2_num_sws_episodes(j) = size(W2_sws_runs{j},1);
	W3_num_sws_episodes(j) = size(W3_sws_runs{j},1);
	W4_num_sws_episodes(j) = size(W4_sws_runs{j},1);

	% Number of REMS episodes
	B_num_rems_episodes(j)  = size(B_rems_runs{j},1);
	W1_num_rems_episodes(j) = size(W1_rems_runs{j},1);
	W2_num_rems_episodes(j) = size(W2_rems_runs{j},1);
	W3_num_rems_episodes(j) = size(W3_rems_runs{j},1);
	W4_num_rems_episodes(j) = size(W4_rems_runs{j},1);

	
% Baseline
	for i=1:size(B_wake_runs{j},1)
    	if B_wake_runs{j}(i,2)-B_wake_runs{j}(i,1)+1 >= 2 
    		B_w_episode_length{j}(i) = B_wake_runs{j}(i,2)-B_wake_runs{j}(i,1)+1;
		end 
	end

	for i=1:size(B_sws_runs{j},1)
    	B_s_episode_length{j}(i) = B_sws_runs{j}(i,2)-B_sws_runs{j}(i,1)+1;
    
	end

	for i=1:size(B_rems_runs{j},1)
		B_r_episode_length{j}(i) = B_rems_runs{j}(i,2)-B_rems_runs{j}(i,1)+1;
	end
	try B_mean_SWS_length(j)  = mean(B_s_episode_length{j}) * 10; catch B_mean_SWS_length(j)  = NaN; end 
	try B_mean_wake_length(j) = mean(B_w_episode_length{j}) * 10; catch B_mean_wake_length(j) = NaN; end
	try B_mean_REM_length(j)  = mean(B_r_episode_length{j}) * 10; catch B_mean_REM_length(j)  = NaN; end 
% W1
	for i=1:size(W1_wake_runs{j},1)
    	if W1_wake_runs{j}(i,2)-W1_wake_runs{j}(i,1)+1 >=2
    		W1_w_episode_length{j}(i) = W1_wake_runs{j}(i,2)-W1_wake_runs{j}(i,1)+1;
    	end 
	end

	for i=1:size(W1_sws_runs{j},1)
    	W1_s_episode_length{j}(i) = W1_sws_runs{j}(i,2)-W1_sws_runs{j}(i,1)+1;
    
	end

	for i=1:size(W1_rems_runs{j},1)
		W1_r_episode_length{j}(i) = W1_rems_runs{j}(i,2)-W1_rems_runs{j}(i,1)+1;
	end
	try W1_mean_SWS_length(j)  = mean(W1_s_episode_length{j}) * 10; catch W1_mean_SWS_length(j)  = NaN; end
	try W1_mean_wake_length(j) = mean(W1_w_episode_length{j}) * 10; catch W1_mean_wake_length(j) = NaN; end 
	try W1_mean_REM_length(j)  = mean(W1_r_episode_length{j}) * 10; catch W1_mean_REM_length(j)  = NaN; end 

% W2
	for i=1:size(W2_wake_runs{j},1)
    	if W2_wake_runs{j}(i,2)-W2_wake_runs{j}(i,1)+1 >=2
    		W2_w_episode_length{j}(i) = W2_wake_runs{j}(i,2)-W2_wake_runs{j}(i,1)+1;
    	end 
	end

	for i=1:size(W2_sws_runs{j},1)
    	W2_s_episode_length{j}(i) = W2_sws_runs{j}(i,2)-W2_sws_runs{j}(i,1)+1;
    
	end

	for i=1:size(W2_rems_runs{j},1)
		W2_r_episode_length{j}(i) = W2_rems_runs{j}(i,2)-W2_rems_runs{j}(i,1)+1;
	end
	try W2_mean_SWS_length(j)  = mean(W2_s_episode_length{j}) * 10; catch W2_mean_SWS_length(j)  = NaN; end
	try W2_mean_wake_length(j) = mean(W2_w_episode_length{j}) * 10; catch W2_mean_wake_length(j) = NaN; end
	try W2_mean_REM_length(j)  = mean(W2_r_episode_length{j}) * 10; catch W2_mean_REM_length(j)  = NaN; end 

% W3
	for i=1:size(W3_wake_runs{j},1)
    	if W3_wake_runs{j}(i,2)-W3_wake_runs{j}(i,1)+1 >=2
    		W3_w_episode_length{j}(i) = W3_wake_runs{j}(i,2)-W3_wake_runs{j}(i,1)+1;
    	end 
	end

	for i=1:size(W3_sws_runs{j},1)
    	W3_s_episode_length{j}(i) = W3_sws_runs{j}(i,2)-W3_sws_runs{j}(i,1)+1;
    
	end

	for i=1:size(W3_rems_runs{j},1)
		W3_r_episode_length{j}(i) = W3_rems_runs{j}(i,2)-W3_rems_runs{j}(i,1)+1;
	end
	try W3_mean_SWS_length(j)  = mean(W3_s_episode_length{j}) * 10; catch W3_mean_SWS_length(j)  = NaN; end 
	try W3_mean_wake_length(j) = mean(W3_w_episode_length{j}) * 10; catch W3_mean_wake_length(j) = NaN; end 
	try W3_mean_REM_length(j)  = mean(W3_r_episode_length{j}) * 10; catch W3_mean_REM_length(j)  = NaN; end 

% W4
	for i=1:size(W4_wake_runs{j},1)
    	if W4_wake_runs{j}(i,2)-W4_wake_runs{j}(i,1)+1 >= 2
    		W4_w_episode_length{j}(i) = W4_wake_runs{j}(i,2)-W4_wake_runs{j}(i,1)+1;
    	end 
	end

	for i=1:size(W4_sws_runs{j},1)
    	W4_s_episode_length{j}(i) = W4_sws_runs{j}(i,2)-W4_sws_runs{j}(i,1)+1;
    
	end

	for i=1:size(W4_rems_runs{j},1)
		W4_r_episode_length{j}(i) = W4_rems_runs{j}(i,2)-W4_rems_runs{j}(i,1)+1;
	end
	try W4_mean_SWS_length(j)  = mean(W4_s_episode_length{j}) * 10; catch W4_mean_SWS_length(j)  = NaN;end 
	try W4_mean_wake_length(j) = mean(W4_w_episode_length{j}) * 10; catch W4_mean_wake_length(j) = NaN;end 
	try W4_mean_REM_length(j)  = mean(W4_r_episode_length{j}) * 10; catch W4_mean_REM_length(j)  = NaN;end 

end  % end of looping over num_simulations



% --- Compute number of transitions between states in each day --- 
for j=1:num_simulations
	W2Sbaseline(j) = length(strfind(state(1:8640,j)','WS')); % baseline
	W2Rbaseline(j) = length(strfind(state(1:8640,j)','WR'));
	S2Wbaseline(j) = length(strfind(state(1:8640,j)','SW'));
	S2Rbaseline(j) = length(strfind(state(1:8640,j)','SR'));
	R2Wbaseline(j) = length(strfind(state(1:8640,j)','RW'));
	R2Sbaseline(j) = length(strfind(state(1:8640,j)','RS'));

	W2S_W1(j) = length(strfind(state(8641:17281,j)','WS')); % W1
	W2R_W1(j) = length(strfind(state(8641:17281,j)','WR'));
	S2W_W1(j) = length(strfind(state(8641:17281,j)','SW'));
	S2R_W1(j) = length(strfind(state(8641:17281,j)','SR'));
	R2W_W1(j) = length(strfind(state(8641:17281,j)','RW'));
	R2S_W1(j) = length(strfind(state(8641:17281,j)','RS'));

	W2S_W2(j) = length(strfind(state(17282:25921,j)','WS')); % W2
	W2R_W2(j) = length(strfind(state(17282:25921,j)','WR'));
	S2W_W2(j) = length(strfind(state(17282:25921,j)','SW'));
	S2R_W2(j) = length(strfind(state(17282:25921,j)','SR'));
	R2W_W2(j) = length(strfind(state(17282:25921,j)','RW'));
	R2S_W2(j) = length(strfind(state(17282:25921,j)','RS'));

	W2S_W3(j) = length(strfind(state(25922:34561,j)','WS')); % W3
	W2R_W3(j) = length(strfind(state(25922:34561,j)','WR'));
	S2W_W3(j) = length(strfind(state(25922:34561,j)','SW'));
	S2R_W3(j) = length(strfind(state(25922:34561,j)','SR'));
	R2W_W3(j) = length(strfind(state(25922:34561,j)','RW'));
	R2S_W3(j) = length(strfind(state(25922:34561,j)','RS'));

	W2S_W4(j) = length(strfind(state(34562:43201,j)','WS')); % W4
	W2R_W4(j) = length(strfind(state(34562:43201,j)','WR'));
	S2W_W4(j) = length(strfind(state(34562:43201,j)','SW'));
	S2R_W4(j) = length(strfind(state(34562:43201,j)','SR'));
	R2W_W4(j) = length(strfind(state(34562:43201,j)','RW'));
	R2S_W4(j) = length(strfind(state(34562:43201,j)','RS'));

end 


% TEST  save raw transition data to .mat file
% if strcmp(shift,'AW')
% 	save('raw_transition_data_AW.mat', 'W2S*', 'W2R*', 'S2W*', 'S2R*', 'R2W*', 'R2S*') 
% elseif strcmp(shift,'RW')
% 	save('raw_transition_data_RW.mat', 'W2S*', 'W2R*', 'S2W*', 'S2R*', 'R2W*', 'R2S*')
% end





% average over the simulations
B_global_av_wake_length   = mean(B_mean_wake_length,'omitnan');
B_global_av_SWS_length    = mean(B_mean_SWS_length,'omitnan');
B_global_av_REM_length    = mean(B_mean_REM_length,'omitnan');
B_global_av_wake_episodes = mean(B_num_wake_episodes,'omitnan');
B_global_av_sws_episodes  = mean(B_num_sws_episodes,'omitnan');
B_global_av_rems_episodes = mean(B_num_rems_episodes,'omitnan');
B_global_av_wake_epochs   = mean(B_num_wake_epochs,'omitnan');
B_global_av_sws_epochs    = mean(B_num_sws_epochs,'omitnan');
B_global_av_rems_epochs   = mean(B_num_rems_epochs,'omitnan');


W1_global_av_wake_length   = mean(W1_mean_wake_length,'omitnan');
W1_global_av_SWS_length    = mean(W1_mean_SWS_length,'omitnan');
W1_global_av_REM_length    = mean(W1_mean_REM_length,'omitnan');
W1_global_av_wake_episodes = mean(W1_num_wake_episodes,'omitnan');
W1_global_av_sws_episodes  = mean(W1_num_sws_episodes,'omitnan');
W1_global_av_rems_episodes = mean(W1_num_rems_episodes,'omitnan');
W1_global_av_wake_epochs   = mean(W1_num_wake_epochs,'omitnan');
W1_global_av_sws_epochs    = mean(W1_num_sws_epochs,'omitnan');
W1_global_av_rems_epochs   = mean(W1_num_rems_epochs,'omitnan');

W2_global_av_wake_length   = mean(W2_mean_wake_length,'omitnan');
W2_global_av_SWS_length    = mean(W2_mean_SWS_length,'omitnan');
W2_global_av_REM_length    = mean(W2_mean_REM_length,'omitnan');
W2_global_av_wake_episodes = mean(W2_num_wake_episodes,'omitnan');
W2_global_av_sws_episodes  = mean(W2_num_sws_episodes,'omitnan');
W2_global_av_rems_episodes = mean(W2_num_rems_episodes,'omitnan');
W2_global_av_wake_epochs   = mean(W2_num_wake_epochs,'omitnan');
W2_global_av_sws_epochs    = mean(W2_num_sws_epochs,'omitnan');
W2_global_av_rems_epochs   = mean(W2_num_rems_epochs,'omitnan');

W3_global_av_wake_length   = mean(W3_mean_wake_length,'omitnan');
W3_global_av_SWS_length    = mean(W3_mean_SWS_length,'omitnan');
W3_global_av_REM_length    = mean(W3_mean_REM_length,'omitnan');
W3_global_av_wake_episodes = mean(W3_num_wake_episodes,'omitnan');
W3_global_av_sws_episodes  = mean(W3_num_sws_episodes,'omitnan');
W3_global_av_rems_episodes = mean(W3_num_rems_episodes,'omitnan');
W3_global_av_wake_epochs   = mean(W3_num_wake_epochs,'omitnan');
W3_global_av_sws_epochs    = mean(W3_num_sws_epochs,'omitnan');
W3_global_av_rems_epochs   = mean(W3_num_rems_epochs,'omitnan');

W4_global_av_wake_length   = mean(W4_mean_wake_length,'omitnan');
W4_global_av_SWS_length    = mean(W4_mean_SWS_length,'omitnan');
W4_global_av_REM_length    = mean(W4_mean_REM_length,'omitnan');
W4_global_av_wake_episodes = mean(W4_num_wake_episodes,'omitnan');
W4_global_av_sws_episodes  = mean(W4_num_sws_episodes,'omitnan');
W4_global_av_rems_episodes = mean(W4_num_rems_episodes,'omitnan');
W4_global_av_wake_epochs   = mean(W4_num_wake_epochs,'omitnan');
W4_global_av_sws_epochs    = mean(W4_num_sws_epochs,'omitnan');
W4_global_av_rems_epochs   = mean(W4_num_rems_epochs,'omitnan');




% std's over the simulations
B_global_std_wake_length   = std(B_mean_wake_length,'omitnan');
B_global_std_SWS_length    = std(B_mean_SWS_length,'omitnan');
B_global_std_REM_length    = std(B_mean_REM_length,'omitnan');
B_global_std_wake_episodes = std(B_num_wake_episodes,'omitnan');
B_global_std_sws_episodes  = std(B_num_sws_episodes,'omitnan');
B_global_std_rems_episodes = std(B_num_rems_episodes,'omitnan');
B_global_std_wake_epochs   = std(B_num_wake_epochs,'omitnan');
B_global_std_sws_epochs    = std(B_num_sws_epochs,'omitnan');
B_global_std_rems_epochs   = std(B_num_rems_epochs,'omitnan');

W1_global_std_wake_length   = std(W1_mean_wake_length,'omitnan');
W1_global_std_SWS_length    = std(W1_mean_SWS_length,'omitnan');
W1_global_std_REM_length    = std(W1_mean_REM_length,'omitnan');
W1_global_std_wake_episodes = std(W1_num_wake_episodes,'omitnan');
W1_global_std_sws_episodes  = std(W1_num_sws_episodes,'omitnan');
W1_global_std_rems_episodes = std(W1_num_rems_episodes,'omitnan');
W1_global_std_wake_epochs   = std(W1_num_wake_epochs,'omitnan');
W1_global_std_sws_epochs    = std(W1_num_sws_epochs,'omitnan');
W1_global_std_rems_epochs   = std(W1_num_rems_epochs,'omitnan');

W2_global_std_wake_length   = std(W2_mean_wake_length,'omitnan');
W2_global_std_SWS_length    = std(W2_mean_SWS_length,'omitnan');
W2_global_std_REM_length    = std(W2_mean_REM_length,'omitnan');
W2_global_std_wake_episodes = std(W2_num_wake_episodes,'omitnan');
W2_global_std_sws_episodes  = std(W2_num_sws_episodes,'omitnan');
W2_global_std_rems_episodes = std(W2_num_rems_episodes,'omitnan');
W2_global_std_wake_epochs   = std(W2_num_wake_epochs,'omitnan');
W2_global_std_sws_epochs    = std(W2_num_sws_epochs,'omitnan');
W2_global_std_rems_epochs   = std(W2_num_rems_epochs,'omitnan');

W3_global_std_wake_length   = std(W3_mean_wake_length,'omitnan');
W3_global_std_SWS_length    = std(W3_mean_SWS_length,'omitnan');
W3_global_std_REM_length    = std(W3_mean_REM_length,'omitnan');
W3_global_std_wake_episodes = std(W3_num_wake_episodes,'omitnan');
W3_global_std_sws_episodes  = std(W3_num_sws_episodes,'omitnan');
W3_global_std_rems_episodes = std(W3_num_rems_episodes,'omitnan');
W3_global_std_wake_epochs   = std(W3_num_wake_epochs,'omitnan');
W3_global_std_sws_epochs    = std(W3_num_sws_epochs,'omitnan');
W3_global_std_rems_epochs   = std(W3_num_rems_epochs,'omitnan');

W4_global_std_wake_length   = std(W4_mean_wake_length,'omitnan');
W4_global_std_SWS_length    = std(W4_mean_SWS_length,'omitnan');
W4_global_std_REM_length    = std(W4_mean_REM_length,'omitnan');
W4_global_std_wake_episodes = std(W4_num_wake_episodes,'omitnan');
W4_global_std_sws_episodes  = std(W4_num_sws_episodes,'omitnan');
W4_global_std_rems_episodes = std(W4_num_rems_episodes,'omitnan');
W4_global_std_wake_epochs   = std(W4_num_wake_epochs,'omitnan');
W4_global_std_sws_epochs    = std(W4_num_sws_epochs,'omitnan');
W4_global_std_rems_epochs   = std(W4_num_rems_epochs,'omitnan');






% Now put the averages and std's together into vectors to make them easy to plot
% Each vector contains 5 entries: B, W1, W2, W3, W4
avg_wake_episodes_vs_time = [B_global_av_wake_episodes W1_global_av_wake_episodes W2_global_av_wake_episodes W3_global_av_wake_episodes W4_global_av_wake_episodes];
avg_SWS_episodes_vs_time  = [B_global_av_sws_episodes W1_global_av_sws_episodes W2_global_av_sws_episodes W3_global_av_sws_episodes W4_global_av_sws_episodes];
avg_REMS_episodes_vs_time = [B_global_av_rems_episodes W1_global_av_rems_episodes W2_global_av_rems_episodes W3_global_av_rems_episodes W4_global_av_rems_episodes];

avg_wake_episode_duration_vs_time = [B_global_av_wake_length    W1_global_av_wake_length    W2_global_av_wake_length    W3_global_av_wake_length    W4_global_av_wake_length];
avg_SWS_episode_duration_vs_time  = [B_global_av_SWS_length    W1_global_av_SWS_length    W2_global_av_SWS_length    W3_global_av_SWS_length    W4_global_av_SWS_length];
avg_REMS_episode_duration_vs_time = [B_global_av_REM_length    W1_global_av_REM_length    W2_global_av_REM_length    W3_global_av_REM_length    W4_global_av_REM_length];

avg_wake_epochs_vs_time = [B_global_av_wake_epochs W1_global_av_wake_epochs W2_global_av_wake_epochs W3_global_av_wake_epochs W4_global_av_wake_epochs];
avg_SWS_epochs_vs_time  = [B_global_av_sws_epochs  W1_global_av_sws_epochs  W2_global_av_sws_epochs  W3_global_av_sws_epochs  W4_global_av_sws_epochs];
avg_REMS_epochs_vs_time = [B_global_av_rems_epochs W1_global_av_rems_epochs W2_global_av_rems_epochs W3_global_av_rems_epochs W4_global_av_rems_epochs];

avg_W2S_transitions = [mean(W2Sbaseline) mean(W2S_W1) mean(W2S_W2) mean(W2S_W3) mean(W2S_W4)];
avg_W2R_transitions = [mean(W2Rbaseline) mean(W2R_W1) mean(W2R_W2) mean(W2R_W3) mean(W2R_W4)];
avg_S2W_transitions = [mean(S2Wbaseline) mean(S2W_W1) mean(S2W_W2) mean(S2W_W3) mean(S2W_W4)];
avg_S2R_transitions = [mean(S2Rbaseline) mean(S2R_W1) mean(S2R_W2) mean(S2R_W3) mean(S2R_W4)];
avg_R2W_transitions = [mean(R2Wbaseline) mean(R2W_W1) mean(R2W_W2) mean(R2W_W3) mean(R2W_W4)];
avg_R2S_transitions = [mean(R2Sbaseline) mean(R2S_W1) mean(R2S_W2) mean(R2S_W3) mean(R2S_W4)];


std_wake_episodes_vs_time = [B_global_std_wake_episodes W1_global_std_wake_episodes W2_global_std_wake_episodes W3_global_std_wake_episodes W4_global_std_wake_episodes];
std_SWS_episodes_vs_time  = [B_global_std_sws_episodes  W1_global_std_sws_episodes  W2_global_std_sws_episodes  W3_global_std_sws_episodes  W4_global_std_sws_episodes];
std_REMS_episodes_vs_time = [B_global_std_rems_episodes W1_global_std_rems_episodes W2_global_std_rems_episodes W3_global_std_rems_episodes W4_global_std_rems_episodes];

std_wake_episode_duration_vs_time = [B_global_std_wake_length   W1_global_std_wake_length   W2_global_std_wake_length   W3_global_std_wake_length   W4_global_std_wake_length];
std_SWS_episode_duration_vs_time  = [B_global_std_SWS_length    W1_global_std_SWS_length    W2_global_std_SWS_length    W3_global_std_SWS_length    W4_global_std_SWS_length];
std_REMS_episode_duration_vs_time = [B_global_std_REM_length    W1_global_std_REM_length    W2_global_std_REM_length    W3_global_std_REM_length    W4_global_std_REM_length];

std_wake_epochs_vs_time = [B_global_std_wake_epochs W1_global_std_wake_epochs W2_global_std_wake_epochs W3_global_std_wake_epochs W4_global_std_wake_epochs];
std_SWS_epochs_vs_time  = [B_global_std_sws_epochs  W1_global_std_sws_epochs  W2_global_std_sws_epochs  W3_global_std_sws_epochs  W4_global_std_sws_epochs];
std_REMS_epochs_vs_time = [B_global_std_rems_epochs W1_global_std_rems_epochs W2_global_std_rems_epochs W3_global_std_rems_epochs W4_global_std_rems_epochs];

std_W2S_transitions = [std(W2Sbaseline) std(W2S_W1) std(W2S_W2) std(W2S_W3) std(W2S_W4)];
std_W2R_transitions = [std(W2Rbaseline) std(W2R_W1) std(W2R_W2) std(W2R_W3) std(W2R_W4)];
std_S2W_transitions = [std(S2Wbaseline) std(S2W_W1) std(S2W_W2) std(S2W_W3) std(S2W_W4)];
std_S2R_transitions = [std(S2Rbaseline) std(S2R_W1) std(S2R_W2) std(S2R_W3) std(S2R_W4)];
std_R2W_transitions = [std(R2Wbaseline) std(R2W_W1) std(R2W_W2) std(R2W_W3) std(R2W_W4)];
std_R2S_transitions = [std(R2Sbaseline) std(R2S_W1) std(R2S_W2) std(R2S_W3) std(R2S_W4)];



% put each into the structs that get returned
averages.num_w_episodes        = avg_wake_episodes_vs_time;
averages.num_sws_episodes      = avg_SWS_episodes_vs_time;  
averages.num_REMS_episodes     = avg_REMS_episodes_vs_time;  
averages.wake_episode_duration = avg_wake_episode_duration_vs_time; 
averages.SWS_episode_duration  = avg_SWS_episode_duration_vs_time;
averages.REMS_episode_duration = avg_REMS_episode_duration_vs_time;
averages.time_in_w    		   = avg_wake_epochs_vs_time*10/60;  % units are minutes
averages.time_in_sws 		   = avg_SWS_epochs_vs_time*10/60;
averages.time_in_rems 		   = avg_REMS_epochs_vs_time*10/60;
averages.W2S_transitions       = avg_W2S_transitions;
averages.W2R_transitions       = avg_W2R_transitions;
averages.S2W_transitions       = avg_S2W_transitions;
averages.S2R_transitions       = avg_S2R_transitions;
averages.R2W_transitions       = avg_R2W_transitions;
averages.R2S_transitions       = avg_R2S_transitions;



stndev.num_w_episodes        = std_wake_episodes_vs_time;
stndev.num_sws_episodes      = std_SWS_episodes_vs_time;  
stndev.num_REMS_episodes     = std_REMS_episodes_vs_time;  
stndev.wake_episode_duration = std_wake_episode_duration_vs_time;  
stndev.SWS_episode_duration  = std_SWS_episode_duration_vs_time;
stndev.REMS_episode_duration = std_REMS_episode_duration_vs_time;
stndev.time_in_w    		 = std_wake_epochs_vs_time*10/60;  % units are minutes
stndev.time_in_sws 		     = std_SWS_epochs_vs_time*10/60;
stndev.time_in_rems 		 = std_REMS_epochs_vs_time*10/60;
stndev.W2S_transitions       = std_W2S_transitions;
stndev.W2R_transitions       = std_W2R_transitions;
stndev.S2W_transitions       = std_S2W_transitions;
stndev.S2R_transitions       = std_S2R_transitions;
stndev.R2W_transitions       = std_R2W_transitions;
stndev.R2S_transitions       = std_R2S_transitions;

if makeplots 
% Finally, plot it all
figure
%subplot(3,1,1)
h1=plot([1:5],avg_wake_episodes_vs_time,'k',[1:5],avg_wake_episodes_vs_time,'ko','MarkerSize',12);
set(h1(2),'MarkerEdgeColor','none','MarkerFaceColor','k')
hold on 
h2=errorbar([1:5],avg_wake_episodes_vs_time,std_wake_episodes_vs_time,'k');
axis([0 5.5 80 240])
ax=gca;
ax.XTick = [1 2 3 4 5];
ax.XTickLabel = {'B', 'W1', 'W2', 'W3','W4'};
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gca,'box','off')
ylabel('Number of wake episodes')
title(shift)

figure
%subplot(3,1,2)
h3=plot([1:5],avg_SWS_episode_duration_vs_time,'k',[1:5],avg_SWS_episode_duration_vs_time,'ko','MarkerSize',12);
set(h3(2),'MarkerEdgeColor','none','MarkerFaceColor','k')
hold on 
h4=errorbar([1:5],avg_SWS_episode_duration_vs_time,std_SWS_episode_duration_vs_time,'k');
axis([0 5.5 100 300])
ax=gca;
ax.XTick = [1 2 3 4 5];
ax.XTickLabel = {'B', 'W1', 'W2', 'W3','W4'};
ax.YTick = [100 150 200 250 300];
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gca,'box','off')
ylabel({'SWS'; '(episode duration, seconds)'})
title(shift)

figure
%subplot(3,1,3)
h5=plot([1:5],avg_REMS_episode_duration_vs_time,'k',[1:5],avg_REMS_episode_duration_vs_time,'ko','MarkerSize',12);
set(h5(2),'MarkerEdgeColor','none','MarkerFaceColor','k')
hold on 
h6=errorbar([1:5],avg_REMS_episode_duration_vs_time,std_REMS_episode_duration_vs_time,'k');
axis([0 5.5 80 140])
ax=gca;
ax.XTick = [1 2 3 4 5];
ax.XTickLabel = {'B', 'W1', 'W2', 'W3','W4'};
ax.YTick = [80 90 100 110 120 130 140];
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gca,'box','off')
ylabel({'REM sleep'; '(episode duration, seconds)'})
title(shift)

end 




