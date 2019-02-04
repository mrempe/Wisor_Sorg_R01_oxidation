function run_costa_model_and_save_output(state)
%  USAGE:  run_costa_model_and_save_output(state)
% main function that solves the SDEs of the Costa et al 2016 model using SDE solver due to Chang 1987 (Mathematics of Computation)
% 
% The simulation is run for 10 minutes currently
% 
% INPUT: state:  either 'wake', 'SWS', or 'REMS'


% first determine if we're using octave or matlab
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
% Octave 
%addpath("Z:/MATLAB/Rempe/Matlab_utilities")
% Matlab
%addpath \\FS1\WisorData\MATLAB\Rempe\Matlab_utilities

% First set up the time vector
dt = 0.1;  % ms  was 0.1
%t=0:dt:86400000;  %86400000 is 24 hours in ms
t=0:dt:600000;  % 3600000 is one hour in ms, 600000 is 10 minutes
n=length(t)-1;

% Now allocate space for Y so it doesn't have to resize every time step
Y = zeros(20,length(t));


% noise parameters from the Costa model
gamma_e	= 70E-3;
dphi	= 2;

v = [zeros(13,1); gamma_e^2*dphi;  gamma_e^2*dphi; zeros(5,1)];  % these are the coeffs of the noise terms, sometimes called g


% Initial Conditions
if strcmp(state,'wake') | strcmp(state,'Wake')
	Y(:,1) = [4;    ...  % F_W was 6
			  1e-3; ...  % F_N
			  1e-3; ...  % F_R
			  0.8;  ...  % C_E
			  1e-3; ...  % C_G
			  1e-3; ...  % C_A
			  0.5;  ...  % h
			  -43;  ...  % V_p was -66
			  -64;  ...  % V_i
			  0;    ...  % s_ep
			  0;    ...  % s_ei
			  0;    ...  % s_gp
			  0;    ...  % s_gi
			  0;    ...  % x_ep (first derivative of s_ep)
			  0;    ...  % x_ei (first derivative of s_ei)
			  0;    ...  % x_gp (first derivative of s_gp)
			  0;    ...  % x_gi (first derivative of s_gi)
			  9.5;  ...  % [Na]
			  1.33; ...  % g_KNa
			  7];        % sigma_p
elseif strcmp(state,'SWS')
	Y(:,1) = [4.7;    ...  % F_W
			  2; ...  % F_N
			  1e-3; ...  % F_R
			  0.6;  ...  % C_E
			  0.6;  ...  % C_G
			  1e-3; ...  % C_A
			  0.9;  ...  % h
			  -45;  ...  % V_p
			  -50;  ...  % V_i
			  0;    ...  % s_ep
			  0;    ...  % s_ei
			  0;    ...  % s_gp
			  0;    ...  % s_gi
			  0;    ...  % x_ep (first derivative of s_ep)
			  0;    ...  % x_ei (first derivative of s_ei)
			  0;    ...  % x_gp (first derivative of s_gp)
			  0;    ...  % x_gi (first derivative of s_gi)
			  28;   ...  % [Na]
			  2;  ...  % g_KNa
			  6.5];        % sigma_p
elseif strcmp(state,'REMS') | strcmp(state,'REM')
	Y(:,1) = [0.93;    ...  % F_W
			  5; ...  % F_N
			  4.9; ...  % F_R
			  0.18;  ...  % C_E
			  0.848; ...  % C_G
			  0.985; ...  % C_A
			  0.705;  ...  % h
			  -55;  ...  % V_p
			  -53;  ...  % V_i
			  0;    ...  % s_ep
			  0;    ...  % s_ei
			  0;    ...  % s_gp
			  0;    ...  % s_gi
			  0;    ...  % x_ep (first derivative of s_ep)
			  0;    ...  % x_ei (first derivative of s_ei)
			  0;    ...  % x_gp (first derivative of s_gp)
			  0;    ...  % x_gi (first derivative of s_gi)
			  9.5;  ...  % [Na]
			  0.2; ...  % g_KNa
			  6.26];        % sigma_p
end



tic
% Start Chang method
for i=1:length(t)-1
	r = randn(2,1);
	alpha  = r(1);
	beta   = 0.5*r(1)+(sqrt(3)/6)*r(2);  
	deltaW = sqrt(dt)*alpha; 
	

	Q        = Y(:,i) + 0.5*dt*F_Costa(Y(:,i));
	Qs       = Y(:,i) + 0.5*dt*F_Costa(Y(:,i)) + (3/2) * v * sqrt(dt)*beta;
	Y(:,i+1) = Y(:,i) + v * deltaW + (1/3)*dt*( F_Costa(Q) + 2*F_Costa(Qs) );
	

	% monitor and report progress
	if i==floor(n/10)
		tenth_time = toc;
		disp(['10% finished... Should be finished at ', char(datetime('now')+9*seconds(tenth_time))])
	end

end


% Now save the results of the simulation in a .mat file that 
% is labeled using the time and day the simulation was run,
% and the sleep state.  
filename = strcat('Model_Output_',date,'__',state,'.mat')
save(filename,'-v7.3')



% TESTING: make a quick plot of V_p
figure
plot(t,Y(8,:))
