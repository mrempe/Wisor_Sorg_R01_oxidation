% main function that calls ode45 on Costa.m using SDE solver due to Chang 1987 (Mathematics of Computation)
% 

% first determine if we're using octave or matlab
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
% Octave 
%addpath("Z:/MATLAB/Rempe/Matlab_utilities")
% Matlab
addpath \\FS1\WisorData\MATLAB\Rempe\Matlab_utilities

clear Y t raw* epoch 

%progressbar

%options = odeset('RelTol',1e-8,'AbsTol',1e-8);
%[t,y] = ode45(@CostaODE,[0 43200000],[6;1e-3;1e-3;0.9;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7]); % use [ 0 86400000] for 24 hours
tic
gamma_e	= 70E-3;
dphi	= 2;

v = [zeros(13,1); gamma_e^2*dphi;  gamma_e^2*dphi; zeros(5,1)];  % these are the coeffs of the noise terms, sometimes called g

% Initial Conditions
Y(:,1) = [6;1e-3;1e-3;0.8;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7];

dt = 0.3;  % ms
%t=0:dt:86400000;  %86400000 is 24 hours in ms
t=0:dt:60000;  % 3600000 is one hour in ms
n=length(t)-1;

% Start Chang method
for i=1:length(t)-1
	r = randn(2,1);
	alpha  = r(1);
	beta   = 0.5*r(1)+(sqrt(3)/6)*r(2);  
	deltaW = sqrt(dt)*alpha; 

	Q        = Y(:,i) + 0.5*dt*F_Costa(Y(:,i));
	Qs       = Y(:,i) + 0.5*dt*F_Costa(Y(:,i)) + (3/2) * v * sqrt(dt)*beta;
	Y(:,i+1) = Y(:,i) + v * deltaW + (1/3)*dt*( F_Costa(Q) + 2*F_Costa(Qs) );
	%progressbar(i/n)
end




toc


% now compute the power in various frequency bands
epoch_length = 10;  % in seconds
freq    	 = 1000/dt;  %10000;  % in Hz
N            = fix(freq*epoch_length);  % number of data points in each epoch
raw_EEG = Y(8,:);
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
% now plot results
figure
plot(t,Y(4,:),'g',t,Y(5,:),'r',t,Y(6,:),'b',t,Y(7,:),'y')  % C_E, C_G, C_A, and h 

figure
plot(t,Y(8,:))   % V_p, a proxy for the EEG signal.  
ax=gca;
ax.YTick = [-90 -60 -30];

figure
plot(epoch.delta)
title('delta power')