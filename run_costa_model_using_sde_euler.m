% main function that calls ode45 on Costa.m using SDE solver due to Chang 1987 (Mathematics of Computation)
%

% first determine if we're using octave or matlab
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% if isOctave 
%     addpath("//FS1/WisorData/MATLAB/Rempe/Matlab_utilities")
% else 
% 	addpath \\FS1\WisorData\MATLAB\Rempe\Matlab_utilities
% end

% MATLAB
addpath \\FS1\WisorData\MATLAB\Rempe\Matlab_utilities
% Octave (next two lines)
%addpath("//FS1/WisorData/MATLAB/Rempe/Matlab_utilities/")
%addpath("//FS1/WisorData/MATLAB/Rempe/Matlab_utilities/SDETools-master/SDETools/")



tic
%options = odeset('RelTol',1e-8,'AbsTol',1e-8);
%[t,y] = ode45(@CostaODE,[0 43200000],[6;1e-3;1e-3;0.9;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7]); % use [ 0 86400000] for 24 hours

gamma_e	= 70E-3;
dphi	= 2;

%v = [zeros(13,1); gamma_e^2*dphi;  gamma_e^2*dphi; zeros(5,1)];

g = @(t,y)[zeros(13,1); gamma_e^2*dphi;  gamma_e^2*dphi; zeros(5,1)];

% Initial Conditions
Y0 = [6; 1e-3; 1e-3; 0.8; 1e-3; 1e-3; 0.5; -66; -64; 0; 0; 0; 0; 0; 0; 0; 0; 9.5; 1.33; 7];

dt = 0.1;  % ms



% --- run simulation in chunks to avoid out-of-memory issues ---
%t=0:dt:86400000;  %86400000 is 24 hours in ms
t1 = 0:dt:1800000;  % 3600000 is one hour in ms
%n = length(t1)-1;
%opts = sdeset('RandSeed',0.5);
% first chunk
Y1 = sde_euler(@(t,y)CostaODE(t,y),@(t,y)g(t,y),t1,Y0);
disp('Finished first run.  Now writing data to disk...')
save('Costa_model_output_chunked.mat','Y1','-v7.3')
Y0 = Y1(end,:);     % make the last step the initial condition for next chunk
clear Y1;			% free up memory
disp('Done writing first chunk to disk. Starting second run...')

% second chunk
t2 = 1800001:dt:3600000;
Y2 = sde_euler(@(t,y)CostaODE(t,y),@(t,y)g(t,y),t2,Y0);
disp('Finished second run. Now writing data to disk...')
save('Costa_model_output_chunked.mat','Y2','-append','-v7.3')
Y0 = Y2(end,:);
clear Y2;






% Now load and concatenate the individual chunks so I can plot them.  

% Y = Y';  % flip so I can reuse plotting code below

% toc
% % now plot results
% figure
% plot(t,Y(4,:),'g',t,Y(5,:),'r',t,Y(6,:),'b',t,Y(7,:),'y')  % C_E, C_G, C_A, and h 

% figure
% plot(t,Y(8,:))   % V_p, a proxy for the EEG signal.  
% ax=gca;
% ax.YTick = [-90 -60 -30];

% V = Y(8,:);
% clear Y 

%save('costa_model_output_sde_euler.mat','t','V')

