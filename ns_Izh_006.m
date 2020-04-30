% ns_Izh_006.m
%24 nov 2016


close all
clear all
clc

% Default values --------------------------------------------------------
%  C = 100; vr = -60; vt = -40; k = 0.7; % parameters used for RS
%  a = 0.03; b = -2; c = -50; d = 100; % neocortical pyramidal neurons
%  vPeak = 35;        % spike cutoff  [mV]
%  dt = 1;            % time step  [ms]

% INPUTS ----------------------------------------------------------------
flagI = 1;           % select input paramters
flagS = 2;           % flags = 5 for spike rates
Sstep = 70;          % step height for input stimulus
S1 = 0.1; S2 = 0.9;  % percentages for step input ON and OFF
dt = 0.5;            
NT = 2000;              % number of time steps

switch flagI
   case 1   % default
       C = 100; vr = -60; vt = -40; k = 0.7; 
       a = 0.03; b = -2; c = -50; d = 100; 
       vPeak = 35;        

   case 2   % miscellaneous
       C = 50; vr = -60; vt = -40; k = 0.7; 
       a = 0.01; b = 5; c = -55; d = 150; 
       vPeak = 10;        
    
    case 3   %   Intrinsically Bursting (IB) Neurons  
        C = 150; vr = -75; vt = -45; k = 1.2; 
        a = 0.01; b = 5; c = -56; d = 130; 
        vPeak = 50; 
        
    case 4   %   Chattering (CH) neuronsdefault
       C = 50; vr = -60; vt = -40; k = 1.5; 
       a = 0.03; b = 1; c = -40; d = 150; 
       vPeak = 35;        
    


end  % switch


% SETUP   ----------------------------------------------------------------
  t = 0:dt:(NT-1) * dt;   % time [ms]

   v = vr*ones(1,NT); u = 0*v; % initial values
   S = zeros(1,NT);            % stimulus

   S(round(S1*NT):round(S2*NT)) =  Sstep ;   % step functon  [70]
   %S = 0.411 .* t;                         % ramp  flagS = 5          
   %S(round(0.25*S1*NT):round(0.75*S2*NT)) = 70;


% SOLVE DE: forward Euler Method =========================================


for m = 1:NT-1  
  vT     = v(m)+ (dt/2) * (k*(v(m) - vr)*(v(m) - vt)-u(m) + S(m))/C;
  v(m+1) = vT + (dt/2)  * (k*(v(m) - vr)*(v(m) - vt)-u(m) + S(m))/C;
 % u(m+1) = u(m) + dt * a*(b*(v(m)-vr)-u(m));
   u(m+1) = u(m) + dt * a*(b*(v(m+1)-vr)-u(m));
    if v(m+1)>= vPeak % a spike is fired!
       v(m)   = vPeak; % padding the spike amplitude
       v(m+1) = c; % membrane voltage reset
       u(m+1) = u(m+1) + d; % recovery variable update
   end; % if
   
end;


% Firing rates: flagF = 1 (yes) flagF = 0 (no) ===========================
   flagF = 1; 
   if flagF == 1
      indexFire = find(v == vPeak);          % indeices for spikes
      indexFire2 = indexFire(2:end);
      indexFire1 = indexFire(1:end-1);
      ISI = dt .* (indexFire2 - indexFire1);  % interspike interval
      fireRate = 1000 ./ISI;                  % fire rate
      f = mean(fireRate);                     % frequency: mean firing rate
      
      F = zeros(1,length(indexFire1));        % average fire rate  [Hz]
      for cc = 2:length(indexFire1)
         F(cc) = (fireRate(cc)+fireRate(cc-1))/2;
      end  % for
      
   disp('Interspike times   ISI   [ms]  ');
   fprintf(' %2.2f ',ISI);
   disp('   ');
   disp('Neuron firing rate  [Hz]   ');
   fprintf(' %2.2f ',fireRate);
   disp('  ');
   fprintf('mean firing rate  f =  %2.2f \n',f);
   disp('  ');
   
   if flagS == 5     % ramp input stimuls
      figure(5)
      fs = 12;
      set(gcf,'units','normalized');
      set(gcf,'position',[0.42 0.05 0.2 0.2]);
      xP = S(indexFire1);  yP = F;
      plot(xP,yP,'b','linewidth',2);
      ylabel('f  [Hz]','fontsize',fs);
      xlabel('S  [pA]','fontsize',fs);
      axis([0 max(xP) 0 1.2*max(yP)])
      grid on
      set(gca,'fontsize',fs);
   end  % if
   
   end  % flagF


% GRAPHICS   =============================================================

figure(1)
   fs = 12;
   set(gcf,'units','normalized');
   set(gcf,'position',[0.02 0.40 0.2 0.25]);
   plot(t, v,'lineWidth',2); 
   xlabel('time  t  [ms]','fontsize',fs);
   ylabel('v  [mV]','fontsize',fs);
   grid on
   set(gca,'fontsize',fs);
   axis([0 1.01*max(t) 1.1*min(v) 1.2*max(v)])

   figure(2)
   fs = 12;
   set(gcf,'units','normalized');
   set(gcf,'position',[0.23, 0.40 0.2 0.25]);
   plot(t, u,'lineWidth',2); 
   xlabel('time  t  [ms]','fontsize',fs);
   ylabel('u  [pA]','fontsize',fs);
   grid on
   set(gca,'fontsize',fs);   
   
figure(3)
   fs = 12;
   set(gcf,'units','normalized');
   set(gcf,'position',[0.44, 0.40 0.2 0.25]);
   plot(v, u,'lineWidth',2); 
   xlabel('v [mV]','fontsize',fs);
   ylabel('u  [pA]','fontsize',fs);
   grid on
   set(gca,'fontsize',fs);     
   
figure(4)
   fs = 12;
   set(gcf,'units','normalized');
   set(gcf,'position',[0.65, 0.40 0.2 0.25]);
   plot(t, S,'lineWidth',2); 
   xlabel('t  [ms]','fontsize',fs);
   ylabel('S  [pA]','fontsize',fs);
   grid on
   set(gca,'fontsize',fs);      
   
   