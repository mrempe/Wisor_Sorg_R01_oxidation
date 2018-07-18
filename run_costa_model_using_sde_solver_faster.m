% run_costa_model_using_sde_solver_faster.m
% 

% first determine if we're using octave or matlab
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
% Octave 
%addpath("Z:/MATLAB/Rempe/Matlab_utilities")
% Matlab
%addpath \\FS1\WisorData\MATLAB\Rempe\Matlab_utilities

clear Y t raw* epoch 
profile on


%options = odeset('RelTol',1e-8,'AbsTol',1e-8);
%[t,y] = ode45(@CostaODE,[0 43200000],[6;1e-3;1e-3;0.9;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7]); % use [ 0 86400000] for 24 hours
tic
gamma_e	= 70E-3;
dphi	= 2;

v = [zeros(13,1); gamma_e^2*dphi;  gamma_e^2*dphi; zeros(5,1)];  % these are the coeffs of the noise terms, sometimes called g


% Initial Conditions
Y = [6;1e-3;1e-3;0.8;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7];



dt = 1;  % ms  was 0.1
%t=0:dt:86400000;  %86400000 is 24 hours in ms
t=0:dt:600000;  % 3600000 is one hour in ms
n=length(t)-1;

% Vp is voltage of pyramidal cells, a proxy for EEG
Vp = zeros(1,length(t));
Vp(1) = Y(8);



tic
% Start Chang method
for i=1:length(t)-1

	r = randn(2,1);
	alpha  = r(1);
	beta   = 0.5*r(1)+(sqrt(3)/6)*r(2);  
	deltaW = sqrt(dt)*alpha; 
	
	% Computing
	F = F_Costa(Y);
	Q        = Y + 0.5*dt*F;
	Qs       = Y + 0.5*dt*F + (3/2) * v * sqrt(dt)*beta;
	Y = Y + v * deltaW + (1/3)*dt*( F_Costa(Q) + 2*F_Costa(Qs) );
	Vp(i+1) = Y(8);
	
	% monitor and report progress
	if i==floor(n/10)
		tenth_time = toc;
		
		disp(['10% finished... Should be finished at ', char(datetime('now')+9*seconds(tenth_time))])
	
	end

end


profile viewer

