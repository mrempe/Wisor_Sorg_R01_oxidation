function [v,avg_mem_potential,network_freq,phi] = PV_model(N,tfinal,dt,gsyn,Iapp_meanIN,Iapp_stdIN)
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


addpath other_models_code   % this is where compute_phi lives



tic
timesteps = tfinal/dt;  
steps_in_one_ms = 1/dt;  
t=0:dt:tfinal;

%re = rand(N,1);             % vertically concatenated. these are column vectors
a  =  0.1*ones(N,1);         % parameter a describes time scale of recovery
b  = -0.1*ones(N,1);      	% parameter b describes sensitivity of u to v
c  =  -70*ones(N,1);         %+15*re.^2;         % parameter c describes after-spike reset value of v
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

% connectivity matrix, connect each cell randomly with all 
% other cells with probability 12%
conn_mat  = rand(N);  % conn_mat: connectivity matrix.  columns are pre-synaptic cells, rows are post-synaptic
%conn_mat = zeros(N);  %TESTING:  remove all connectivity

locs_low  = find(conn_mat<=0.12);
locs_high = find(conn_mat>0.12);
conn_mat(locs_low)  = 1;
conn_mat(locs_high) = 0;


[post_syn_cells,pre_syn_cells] = find(conn_mat);
post_syn_cells = unique(post_syn_cells);
pre_syn_cells  = unique(pre_syn_cells);



for i=post_syn_cells'
	pre_syn_connections{i} = find(conn_mat(i,:));  % cell array where each entry corresponds to a post-synaptic cell. Each entry lists pre-synaptic cells, if any, that connect to this cell
end


v = 10*rand(N,1)-65*ones(N,1); %-65*ones(N,1);    % Initial values of v, all -65 mV
u = b.*v;            % Initial values of u done in original Izhevich code
%u = zeros(size(v));  % try initializing u to zero
s = zeros(length(v),timesteps);
firings=[];             % spike timings

T    = zeros(N,timesteps);		  % rows correspond to cell numbers, columns to time steps
Isyn = zeros(N,timesteps);			% initialize vector for synaptic input

for step=1:timesteps           
  vplot(step)=v(1,step);
  v_with_spikes(:,step)=v(:,step);    % Before you reset v, save it
  uplot(step)=u(1);
  Iapp = randn(N,1).*Iapp_std + Iapp_mean;
	fired    = find(v(:,step)>=2.5);    			 % indices of cells that spiked
  fired_cells_indx(:,step) = v(:,step)>=2.5;            % logical array of whether cell fire 
  %firings  = [firings; t(step)+0*fired,fired]; % add a row to firings for every cell that firedwith two elements: time fired and cell number
                                                % this works, but is slow.  You can uncomment previous line if needed. 

  if ~isempty(fired)                        % if at least one cell fired this timestep, reset v,u,T, and k for those cells
    v(fired,step) = c(fired);               % after spike reset v
    u(fired)      = u(fired)+d(fired);	          % after spike reset u
    T(fired,step:step+steps_in_one_ms) = 1;
    % v_high   = find(v(fired,step)>=vt);
    % v_low    = find(v(fired,step)<vt);
    % k(v_high)= k_high;
    % k(v_low) = k_low;
    k(v(fired,step)>=vt) = k_high;  %Faster? logical indexing
    k(v(fired,step)<vt)  = k_low;  
  else
    k(:) = k_low;                             % if no cells fired this timestep, don't reset v, or u. Just set k to low value
  end
  
  
  for j = post_syn_cells'
    if any(T(pre_syn_connections{j},step))  % if any presyanptic cell has fired in the last ms
      s(j,step+1) = s_inf + (s(j,step)-s_inf)*exp(-dt/tau_s);
    else
      s(j,step+1) = s(j,step)*exp(-beta*dt);
    end
  end



  conn_mat_pre_syn_fired_cols = conn_mat(:,fired);
  [rows,cols]=find(conn_mat_pre_syn_fired_cols);
  if ~isempty(rows)
  	for i=1:length(rows)
		  Isyn(rows(i),step) = gsyn*s(rows(i),step)*(v(rows(i))-esyn);
  	end
  end
  v(:,step+1) = v(:,step)+(dt./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn(:,step));
  u = u+dt.*a.*(b.*(v(:,step)-vr)-u); 
             
  



  splot(step) = s(5,step);
end  % end of looping through timesteps

avg_mem_potential = mean(v_with_spikes,1);


% compute network frequency following Ferguson et al 2013 
% (freq at which there is a spectral peak in the last 500ms of population activity)
fs = (1/dt)*1000;     % sampling freq. in Hz  1000 to get from ms to sec
signal = avg_mem_potential(100001:end)-mean(avg_mem_potential(100001:end));  % remove mean to detrend
[pxx,f] = pwelch(signal,fs/2,round(0.6*fs/2),200000,fs);
max_freq_indx = find(pxx==max(pxx));
network_freq = f(max_freq_indx)

% compute network coherency
phi = compute_phi(v,N,dt)



toc
figure
%plot(firings(:,1),firings(:,2),'.')  %slow but works well
spy(fired_cells_indx)
title(['Raster Plot gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(250))])
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
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(250))])
xlabel('Time (ms)')
ylabel('Synaptic gating variable s'  )

figure
plot(t(1:end-1),avg_mem_potential)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential)')

figure
plot(f,10*log10(pxx),f,10*log10(pxx),'x')
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(250))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz)')
axis([0 350 ylim])



