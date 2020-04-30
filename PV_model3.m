function [v,avg_mem_potential,network_freq,phi,spectral_power] = PV_model3(N,tfinal,dt,gsyn,Iapp_meanIN,Iapp_stdIN)
% 
% Usage: v=Izhikevich_modified(N,tfinal,gsyn,Iapp_mean,Iapp_std)
%
% This function attempts to reproduce simulation results from
% Ferguson et al 2013.  A network of PV-positive cells using 
% a framework from Izhikevich 2003
%
% INPUTS defaults listed in brackets: 	
%     N:  		      # of cells  [500]
% 		tfinal: 	    length of simulation (in ms) [1500]
%     dt:           timestep (in ms) [0.01]
%			gsyn: 		    synapse conductance (in nS) [150 NOTE typo in Ferguson et al 2013]
%			Iapp_meanIN:  mean applied current (in pA) [600]
%			Iapp_stdIN:	  std of applied current (in pA) [12]
%
% OUTPUT:	v: 	voltage of all cells in the network.  
%         avg_mem_potential: membrane potential averaged across cells (contains as many elements as timesteps)
%         network_freq:  frequency at which there is a spectral peak during last 500ms of simulation
%         phi:            average cross-correlation between all spike trains (whole simulation)
tic



addpath other_models_code   % this is where compute_phi lives



timesteps = tfinal/dt;  
steps_in_one_ms = 1/dt;  
t=0:dt:tfinal;
last500ms_indx = find(abs(t-(tfinal-500))<1e-12); % just the one index corresponding to time t=500 ms before end
last_500ms_indices = find(t>=tfinal-500);
last_500ms_indices = last_500ms_indices(1:end-1);

%re = rand(N,1);             % vertically concatenated. these are column vectors
a  =  0.1*ones(N,1);         % parameter a describes time scale of recovery
b  = -0.1*ones(N,1);      	% parameter b describes sensitivity of u to v
c  =  -67*ones(N,1);         %+15*re.^2;         % parameter c describes after-spike reset value of v
d  =   10*ones(N,1);          %0.1*ones(N,1);   %8-6*re.^2;           	 % parameter d describes after-spike reset value of u
%S  =  0.5*rand(N,N);  
Cm =   90*ones(N,1);
k  =  1.7*ones(N,1);
vr = -60.6*ones(N,1);
vt = -43.1;  %*ones(N,1);
alpha = (1/0.27);
beta  = (1/1.8);
s_inf = alpha/(alpha+beta);
tau_s = 1/(alpha+beta);
esyn  = -85;   %*ones(N,1);   
k_low = 1.7;
k_high = 14;


Iapp_std  = Iapp_stdIN*ones(N,1);
Iapp_mean = Iapp_meanIN*ones(N,1);
Iapp  = randn(N,1).*Iapp_std + Iapp_mean;

% connectivity matrix, connect each cell randomly with all 
% other cells with probability 12%
conn_mat  = rand(N);  % conn_mat: connectivity matrix.  columns are pre-synaptic cells, rows are post-synaptic
%conn_mat = zeros(N);  %TESTING:  remove all connectivity

locs_low  = find(conn_mat<=0.12);
locs_high = find(conn_mat>0.12);
conn_mat(locs_low)  = 1;
conn_mat(locs_high) = 0;
conn_mat  = conn_mat-diag(diag(conn_mat));  % remove the diagonal entries so cells don't synapse onto themselves
mapping = cell(N,1);          % cell array listing all presynaptic cells for each post-syn cell
for i=1:N
  mapping{i} = find(conn_mat(i,:));
end

% [post_syn_cells,pre_syn_cells] = find(conn_mat);
% post_syn_cells = unique(post_syn_cells);
% pre_syn_cells  = unique(pre_syn_cells);



% for i=post_syn_cells'
% 	pre_syn_connections{i} = find(conn_mat(i,:));  % cell array where each entry corresponds to a post-synaptic cell. Each entry lists pre-synaptic cells, if any, that connect to this cell
% end

% inline functions for RK4 or Adams-Bashforth
% fv = inline('(1/Cm)*(k.*(v-vr).*(v-vt)-u+Iapp-Isyn)','v','Cm','k','vr','vt','u','Isyn','Iapp');
% fu = inline('a.*(b.*(v-vr)-u)','u','v','a','b','vr');

%fv = @(v,Cm,k,vr,vt,u,Isyn,Iapp) (1/Cm)*(k.*(v-vr).*(v-vt)-u+Iapp-Isyn);
%fu = @(u,v,a,b,vr) a.*(b.*(v-vr)-u);


v = 10*rand(N,1)-65*ones(N,1); %-65*ones(N,1);    % Initial values of v, all -65 mV
u = b.*v;            % Initial values of u done in original Izhevich code
%u = zeros(size(v));  % try initializing u to zero
s = zeros(length(v),timesteps);
firings=[];             % spike timings

T    = zeros(N,timesteps);		  % rows correspond to cell numbers, columns to time steps
Isyn = zeros(N,1);			% initialize vector for synaptic input


% Bootstrap: Take one step of FE before starting Adams-Bashforth method 
s_now = s(:,1);
CS = zeros(N,1);
for i=1:N
  CS(i) = sum(s_now(mapping{i}));
end

Isyn = (gsyn*CS).*(v(:,1)-esyn);
uold = u;
Isyn_old = Isyn;
k_old = k;
v(:,2) = v(:,1)+(dt./Cm).*(k.*(v(:,1)-vr).*(v(:,1)-vt)-u+Iapp-Isyn);
unew = u+dt.*a.*(b.*(v(:,1)-vr)-u); 
vplot(1)           = v(1,1);
v_with_spikes(:,1) = v(:,1);    
uplot(1)           = u(1);



% --- COMPUTING LOOP ---------------------------------------------------------------------------------
for step=2:timesteps           
  vplot(step)           = v(1,step);
  v_with_spikes(:,step) = v(:,step);    % Before you reset v, save it
  uplot(step)           = u(1);
  
	fired = find(v(:,step)>=2.5);    			 % indices of cells that spiked
  fired_cells_indx(:,step) = v(:,step)>=2.5;            % logical array of whether cell fire 
  %firings  = [firings; t(step)+0*fired,fired]; % add a row to firings for every cell that firedwith two elements: time fired and cell number
                                                % this works, but is slow.  You can uncomment previous line if needed. 
  v_high   = find(v(:,step)>vt);
  v_low    = find(v(:,step)<=vt);
  k(v_high)= k_high;
  k(v_low) = k_low;


  if ~isempty(fired)                        % if at least one cell fired this timestep, reset v,u,T, and k for those cells
    v(fired,step) = c(fired);               % after spike reset v
    u(fired)      = u(fired)+d(fired);	          % after spike reset u
    T(fired,step:step+steps_in_one_ms) = 1;
  end
  % figure(47)
  % spy(T(:,1:step+steps_in_one_ms))
  % axis([0 step+steps_in_one_ms 0 501])
 
  % OLD WAY: 
  % for j = post_syn_cells'
  %   if any(T(pre_syn_connections{j},step))  % if any presyanptic cell has fired in the last ms
  %     s(j,step+1) = s_inf + (s(j,step)-s_inf)*exp(-dt/tau_s);
  %   else
  %     s(j,step+1) = s(j,step)*exp(-beta*dt);
  %   end
  % end

% OLD WAY: 
  % conn_mat_pre_syn_fired_cols = conn_mat(:,fired);
  % [rows,cols]=find(conn_mat_pre_syn_fired_cols);
  % if ~isempty(rows)
  % 	for i=1:length(rows)
		%   Isyn(rows(i),step) = gsyn*s(rows(i),step)*(v(rows(i))-esyn);
  % 	end
  % end

% NEW WAY: (works great, but is a bit slow)
%Isyn = (gsyn*conn_mat*s(:,step)).*(v(:,step)-esyn);

% Newer (and faster way)
s_now = s(:,step);
CS = zeros(N,1);
for i=1:N
  CS(i) = sum(s_now(mapping{i}));
end

Isyn = (gsyn*CS).*(v(:,step)-esyn);

% Forward Euler (works 3.21.2019)
%v(:,step+1) = v(:,step)+(dt./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn);
%u = u+dt.*a.*(b.*(v(:,step)-vr)-u); 
% FE new way: (works 3.22.2019 1:55PM)
% v(:,step+1) = v(:,step)+dt*fv(v(:,step),Cm,k,vr,vt,u,Isyn,Iapp);
% u = u+dt*fu(u,v(:,step),a,b,vr);


% Second order Adams-Bashforth (working as of 3.22.2019 1:57PM)
% v(:,step+1) = v(:,step) + (dt/2).*(3*fv(v(:,step),Cm,k,vr,vt,u,Isyn,Iapp)-fv(v(:,step-1),Cm,k_old,vr,vt,uold,Isyn_old,Iapp));
% unew        = u         + (dt/2).*(3*fu(u,v(:,step),a,b,vr)-fu(uold,v(:,step-1),a,b,vr));
% uold = u;
% u = unew;
% Isyn_old = Isyn;
% k_old = k;

% RK4 (problem here may be evaluating Isyn and u at the intermediate time steps)
uhalf = u+(dt/2).*a.*(b.*(v(:,step)-vr)-u);
ufull = u+dt.*a.*(b.*(v(:,step)-vr)-u);
fv1 = fv(v(:,step),           Cm,k,vr,vt,u,Isyn,Iapp,0);
fv2 = fv(v(:,step)+(dt/2)*fv1,Cm,k,vr,vt,uhalf,(gsyn*CS).*(v(:,step)+(dt/2)*fv1-esyn),Iapp,0);
fv3 = fv(v(:,step)+(dt/2)*fv2,Cm,k,vr,vt,uhalf,(gsyn*CS).*(v(:,step)+(dt/2)*fv2-esyn),Iapp,0);
fv4 = fv(v(:,step)+(dt)*fv3,  Cm,k,vr,vt,ufull,(gsyn*CS).*(v(:,step)+(dt)*fv3-esyn),Iapp,0);
v(:,step+1) = v(:,step) + (dt/6).*(fv1+2*fv2+2*fv3+fv4);

vhalf = v(:,step)+((dt/2)./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn);
vfull = v(:,step)+(dt./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn);
fu1 = fu(u,           v(:,step),a,b,vr);
fu2 = fu(u+(dt/2)*fu1,vhalf,a,b,vr);
fu3 = fu(u+(dt/2)*fu2,vhalf,a,b,vr);
fu4 = fu(u+(dt)*fu3,  vfull,a,b,vr);
u = u + (dt/6).*(fu1+2*fu2+2*fu3+fu4);


rising_locs  = find(T(:,step));
falling_locs = find(T(:,step)==0);

s(rising_locs,step+1)  =  s_inf + (s(rising_locs,step)-s_inf)*exp(-dt/tau_s);
s(falling_locs,step+1) = s(falling_locs,step)*exp(-beta*dt);
% Old way, too slow: 
% for i=1:N % loop over cells
%   if T(i,step)
%     s(i,step+1) = s_inf + (s(i,step)-s_inf)*exp(-dt/tau_s);
%   else
%     s(i,step+1) = s(i,step)*exp(-beta*dt);
%   end
% end



  splot(step) = s(5,step);

  if step ==round(timesteps/4)
    et = toc;
    disp(['25% finished.... Should be finished at ',num2str(datestr(now + 3*seconds(et)))])
  elseif step == round(timesteps/2)
    disp('50% finished...')
  elseif step == round(0.75*timesteps)
    disp('75% finished...')
  end 




end  % end of looping through timesteps

avg_mem_potential = mean(v_with_spikes,1);


% compute network frequency following Ferguson et al 2013 
% (freq at which there is a spectral peak in the last 500ms of population activity)
fs      = (1/dt)*1000;     % sampling freq. in Hz  1000 to get from ms to sec
signal  = avg_mem_potential(last_500ms_indices)-mean(avg_mem_potential(last_500ms_indices));  % remove mean to detrend
[pxx,f] = pwelch(signal,fs/2,round(0.6*fs/2),2*fs,fs);
max_freq_indx = find(pxx==max(pxx));
network_freq  = f(max_freq_indx)
gammaIdx    = (f>=80 & f<=100);
theta530Idx = (f>=5 & f<=30);
spectral_power = struct;
spectral_power.gamma = sum(pxx(gammaIdx,:),1);
spectral_power.theta = sum(pxx(theta530Idx,:),1); 


% compute network coherency
phi = compute_phi(v(:,last_500ms_indices),N,dt);
phi2 = compute_phi2(fired_cells_indx,N,dt,last500ms_indx,network_freq)
phi = phi2;  


figure
%plot(firings(:,1),firings(:,2),'.')  %slow but works well
spy(fired_cells_indx)
title(['Raster Plot gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(5))])
ylabel('Cell number')
%xlabel('Time (ms)')  %use this if using plot(firings(:,1),firings(:,2),'.')
xlabel('Timesteps')
axis normal


figure
plot(t(1:end-1),v_with_spikes(1,:))
title(['Output of cell 1, gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Voltage (mV)')


figure
plot(t(1:end-1),splot)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(5))])
xlabel('Time (ms)')
ylabel('Synaptic gating variable s'  )

figure
plot(t(1:end-1),avg_mem_potential)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential)')

figure
plot(f,10*log10(pxx),f,10*log10(pxx),'x')
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(5))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz)')
axis([0 350 ylim])




