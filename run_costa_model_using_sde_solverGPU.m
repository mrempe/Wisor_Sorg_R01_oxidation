% main function that calls ode45 on Costa.m using SDE solver due to Chang 1987 (Mathematics of Computation)
% USING THE GPU

% first determine if we're using octave or matlab
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
% Octave 
%addpath("Z:/MATLAB/Rempe/Matlab_utilities")
% Matlab
addpath \\FS1\WisorData\MATLAB\Rempe\Matlab_utilities


progressbar

%options = odeset('RelTol',1e-8,'AbsTol',1e-8);
%[t,y] = ode45(@CostaODE,[0 43200000],[6;1e-3;1e-3;0.9;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7]); % use [ 0 86400000] for 24 hours
tic
gamma_e	= 70E-3;
dphi	= 2;

v = [zeros(13,1); gamma_e^2*dphi;  gamma_e^2*dphi; zeros(5,1)];  % these are the coeffs of the noise terms, sometimes called g

% Initial Conditions
Y(:,1) = [6;1e-3;1e-3;0.8;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7];

Q  = zeros(size(Y));
Qs = zeros(size(Y));

% GPU versions of each variable
GY  = gpuArray(Y);
GQ  = gpuArray(Q);
GQs = gpuArray(Qs);


dt = 1;  % ms
t=0:dt:86400000;  %86400000 is 24 hours in ms
%t=0:dt:60000;  % 3600000 is one hour in ms
n=length(t)-1;

% Start Chang method
for i=1:length(t)-1
	r = randn(2,1);
	alpha  = r(1);
	beta   = 0.5*r(1)+(sqrt(3)/6)*r(2);  
	deltaW = sqrt(dt)*alpha; 

	GQ         = GY(:,i) + 0.5*dt*F_Costa(GY(:,i));
	GQs        = GY(:,i) + 0.5*dt*F_Costa(GY(:,i)) + (3/2) * v * sqrt(dt)*beta;
	GY(:,i+1) = GY(:,i) + v * deltaW + (1/3)*dt*( F_Costa(Q) + 2*F_Costa(Qs) );
	progressbar(i/n)
end