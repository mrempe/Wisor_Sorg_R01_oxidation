function power = compute_power_in_freq_bands(av_mem_pot,dt,last_500ms_indices)
%
% USAGE: gamma_power = compute_gamma_power(av_mem_pot,dt)
%
% INPUTS: 
% av_mem_pot: 	average membrane potential, one of the outputs of 
%				Izhikevich_modified.m
% dt:			time step size (in ms)


% compute the power in various frequency bands
epoch_length = 0.01;  % in seconds
freq    	 = 1000*(1/dt);  %10000;  % in Hz since dt is in ms
f=0.5:0.5:150;
N            = fix(freq*epoch_length);  % number of data points in each epoch
raw_signal = av_mem_pot(last_500ms_indices)-mean(av_mem_pot(last_500ms_indices)); % use only last 0.5 seconds of signal and detrend
raw_signal_truncated = raw_signal(1:N*floor(numel(raw_signal)/N));
raw_signal_reshaped  = reshape(raw_signal_truncated,N,length(raw_signal_truncated)/N);
[spectra_new,freqs_new] = pwelch(raw_signal_reshaped,hamming(floor(N/2)),[],f,freq);
[spectra_new,freqs_new] = pwelch(raw_signal_reshaped,[],[],f,freq);

		
deltaIdx    = (freqs_new>=1  & freqs_new<=4);		%find the indices first
thetaIdx    = (freqs_new>=5  & freqs_new<=9);
gammaIdx    = (freqs_new>=80 & freqs_new<=100);
alphaIdx    = (freqs_new>=8  & freqs_new<=15);
epoch = struct;

epoch.delta = sum(spectra_new(deltaIdx,:),1);   %creates a row vector of delta power in each epoch (may be nice if you want to plot power vs time)
epoch.theta = sum(spectra_new(thetaIdx,:),1);
epoch.gamma = sum(spectra_new(gammaIdx,:),1);
epoch.alpha = sum(spectra_new(alphaIdx,:),1);

power = struct;
power.delta = mean(epoch.delta);
power.theta = mean(epoch.theta);
power.gamma = mean(epoch.gamma);
power.alpha = mean(epoch.alpha); 

disp(['Average gamma power: ', num2str(mean(epoch.gamma))])
