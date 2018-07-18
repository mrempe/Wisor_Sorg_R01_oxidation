% main function that calls ode45 on Costa.m
% 
addpath \\FS1\WisorData\MATLAB\Rempe\Matlab_utilities


%options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,y] = ode45(@CostaODE,[0 43200000],[6;1e-3;1e-3;0.9;1e-3;1e-3;0.5;-66;-64;0;0;0;0;0;0;0;0;9.5;1.33;7]); % use [ 0 86400000] for 24 hours







% now plot results
figure
plot(t,y(:,4),'g',t,y(:,5),'r',t,y(:,6),'b',t,y(:,7),'y')  % C_E, C_G, C_A, and h 

figure
plot(t,y(:,8))   % V_p, a proxy for the EEG signal.  
ax=gca;
ax.YTick = [-90 -60 -30];