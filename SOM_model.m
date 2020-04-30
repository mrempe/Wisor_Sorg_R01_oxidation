function [v,avg_mem_potential,network_freq,phi2] = SOM_model(N,tfinal,dt,gsyn,Iapp_meanIN,Iapp_stdIN)
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



timesteps = tfinal/dt;  
steps_in_one_ms = 1/dt;  
t=0:dt:tfinal;
last500ms_indx = find(abs(t-(tfinal-500))<1e-12);
last_500ms_indices = find(t>=tfinal-500);
last_500ms_indices = last_500ms_indices(1:end-1);  % remove the last one so sizes work


%re = rand(N,1);             % vertically concatenated. these are column vectors
a  =  0.001*ones(N,1);         % parameter a describes time scale of recovery
b  = 0.2*ones(N,1);      	% parameter b describes sensitivity of u to v
c  =  -64*ones(N,1);         %+15*re.^2;         % parameter c describes after-spike reset value of v
d  =   8*ones(N,1);          %0.1*ones(N,1);   %8-6*re.^2;           	 % parameter d describes after-spike reset value of u
%S  =  0.5*rand(N,N);  
Cm =   100*ones(N,1);
k  =  1.7*ones(N,1);
vr = -64*ones(N,1);  % -64 from Hu et al 2011
vt = -35;  %*ones(N,1);  % threshold set to -35 mV from Hu et al 2011
alpha = (1/0.27);
beta  = (1/1.8);
s_inf = alpha/(alpha+beta);
tau_s = 1/(alpha+beta);
esyn  = -85;   %*ones(N,1);   
k_low = 0.34; %1.4;  % my guess for SOM cells
k_high = 1.5; %1.4;


% ----- Gap junction parameters -------------------------------------------------------
g_gap = 2; %2.1 +- 2.3 nS according to Gibson 1999 recorded in rat ;  % in nS  
% -------------------------------------------------------------------------------------




Iapp_std  = Iapp_stdIN*ones(N,1);
Iapp_mean = Iapp_meanIN*ones(N,1);
Iapp  = randn(N,1).*Iapp_std + Iapp_mean;

% connectivity matrix, connect each cell randomly with all 
% other cells with probability 12%
conn_mat  = rand(N);  % conn_mat: connectivity matrix.  columns are pre-synaptic cells, rows are post-synaptic
%conn_mat = zeros(N);  %TESTING:  remove all connectivity

conn_mat = conn_mat-diag(diag(conn_mat));  % remove the diagonal entries so cells don't synapse onto themselves
locs_low  = find(conn_mat<=0.12);
locs_high = find(conn_mat>0.12);
conn_mat(locs_low)  = 1;
conn_mat(locs_high) = 0;

mapping = cell(N,1);          % cell array listing all presynaptic cells for each post-syn cell
for i=1:N
  mapping{i} = find(conn_mat(i,:));  %his just has locations.  Make it include values too somehow,
end
conductance = cell(N,1);
for row=1:N
  for i=1:length(mapping{row})
    conductance{row}(i) = conn_mat(row,mapping{row}(i));
  end
end                               % conductance contains the values of all the gsyn, not just locations.  
cells_receiving_syn_input = find(~cellfun(@isempty,mapping));

% [post_syn_cells,pre_syn_cells] = find(conn_mat);
% post_syn_cells = unique(post_syn_cells);
% pre_syn_cells  = unique(pre_syn_cells);



% for i=post_syn_cells'
% 	pre_syn_connections{i} = find(conn_mat(i,:));  % cell array where each entry corresponds to a post-synaptic cell. Each entry lists pre-synaptic cells, if any, that connect to this cell
% end


% ------------ GAP JUNCTION CONNECTIVITY FOR SOM POPULATION -------------------------------------
gap_conn_mat = rand(N);
locs_low     = find(gap_conn_mat<=0.93);  % 93% connectivity between SOM cells
locs_high    = find(gap_conn_mat> 0.93);
gap_conn_mat(locs_low)  = 1;
gap_conn_mat(locs_high) = 0;
gap_conn_mat = gap_conn_mat-triu(gap_conn_mat)+(tril(gap_conn_mat))';  % make it symmetric
gap_conn_mat = gap_conn_mat.*(1-eye(N));                                % remove main diagonal
%gap_conn_mat = g_gap*gap_conn_mat;
% -----------------------------------------------------------------------------------------------




v = 10*rand(N,1)-64*ones(N,1); %-65*ones(N,1);    % Initial values of v, all -65 mV
u = b.*v;            % Initial values of u done in original Izhevich code
%u = zeros(size(v));  % try initializing u to zero
s = zeros(length(v),timesteps);
firings=[];             % spike timings

T    = zeros(N,timesteps);		  % rows correspond to cell numbers, columns to time steps
Isyn = zeros(N,timesteps);			% initialize vector for synaptic input
Igap = zeros(N,1);


s_now = s(:,1);
CS = zeros(N,1);
for i=1:length(cells_receiving_syn_input) 
  cond_this_cell = conductance{cells_receiving_syn_input(i)};   % this should give a row vector of gsyn values for all cells presynaptic to cell i
  CS(i) = cond_this_cell*(s_now(mapping{cells_receiving_syn_input(i)}));  % dot product: row vector multiplied by a column vec 
end


%Isyn = (gsyn*CS).*(v(:,1)-esyn); % Old way, when gsyn was constant for whole network
Isyn = (CS).*(v(:,1)-esyn);




tic 
disp('Done initializing.  Now updating the network....')
for step=1:timesteps           
  vplot(step)           = v(1,step);
  v_with_spikes(:,step) = v(:,step);    % Before you reset v, save it
  uplot(step)           = u(1);
  
	fired = find(v(:,step)>=11);    			 % indices of cells that spiked. I'm using 11 mV because I measured spike height in Hu et al 2011.  
  fired_cells_indx(:,step) = v(:,step)>=11;            % logical array of whether cell fire 
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
 


  %Isyn(:,step) = (gsyn*conn_mat*s(:,step)).*(v(:,step)-esyn); % works, but a bit slow
  s_now = s(:,step);
  CS = zeros(N,1);
  for i=1:length(cells_receiving_syn_input)
    %CS(i) = sum(s_now(mapping{i}));
    cond_this_cell = conductance{cells_receiving_syn_input(i)};   % this should give a row vector of gsyn values for all cells presynaptic to cell i
    CS(i) = cond_this_cell*(s_now(mapping{cells_receiving_syn_input(i)}));  % dot product: row vector multiplied by a column vec 
  end
  %Isyn(:,step) = (gsyn*CS).*(v(:,step)-esyn);
  Isyn(:,step) = (CS).*(v(:,step)-esyn);

  % --- gap junctions for SOM cells ----------------------
  % for i=1:N
  %   Igap(i) = 0;
  %   for j=1:N
  %     Igap(i) = Igap(i) + gap_conn_mat(i,j)*(v(i,step)-v(j,step));
  %   end
  % end
  vdiffs = gap_conn_mat.*(repmat(v(:,step)',N,1) - repmat(v(:,step),1,N));
  Igap = g_gap*sum(vdiffs,2);



  % Forward Euler (works)
  % v(:,step+1) = v(:,step)+(dt./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn(:,step));
  % u = u+dt.*a.*(b.*(v(:,step)-vr)-u); 

% RK4 
  uhalf = u+(dt/2).*a.*(b.*(v(:,step)-vr)-u);
  ufull = u+dt.*a.*(b.*(v(:,step)-vr)-u);
  fv1 = fv(v(:,step),           Cm,k,vr,vt,u,Isyn(:,step),Iapp,Igap);
  fv2 = fv(v(:,step)+(dt/2)*fv1,Cm,k,vr,vt,uhalf,(CS).*(v(:,step)+(dt/2)*fv1-esyn),Iapp,Igap);
  fv3 = fv(v(:,step)+(dt/2)*fv2,Cm,k,vr,vt,uhalf,(CS).*(v(:,step)+(dt/2)*fv2-esyn),Iapp,Igap);
  fv4 = fv(v(:,step)+(dt)*fv3,  Cm,k,vr,vt,ufull,(CS).*(v(:,step)+(dt)*fv3-esyn),Iapp,Igap);
  v(:,step+1) = v(:,step) + (dt/6).*(fv1+2*fv2+2*fv3+fv4);

  vhalf = v(:,step)+((dt/2)./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn(:,step));
  vfull = v(:,step)+(dt./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn(:,step));
  fu1 = fu(u,       v(:,step),a,b,vr);
  fu2 = fu(u+(dt/2)*fu1,vhalf,a,b,vr);
  fu3 = fu(u+(dt/2)*fu2,vhalf,a,b,vr);
  fu4 = fu(u+(dt)*fu3,  vfull,a,b,vr);
  u = u + (dt/6).*(fu1+2*fu2+2*fu3+fu4);


  rising_locs  = find(T(:,step));
  falling_locs = find(T(:,step)==0);

  s(rising_locs,step+1)  =  s_inf + (s(rising_locs,step)-s_inf)*exp(-dt/tau_s);
  s(falling_locs,step+1) = s(falling_locs,step)*exp(-beta*dt);
   


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
fs = (1/dt)*1000;     % sampling freq. in Hz  1000 to get from ms to sec
signal = avg_mem_potential(last_500ms_indices)-mean(avg_mem_potential(last_500ms_indices));  % remove mean to detrend
[pxx,f] = pwelch(signal,fs/2,round(0.6*fs/2),2*fs,fs);
max_freq_indx = find(pxx==max(pxx));
network_freq = f(max_freq_indx)

% compute network coherency
phi = compute_phi(v(:,last_500ms_indices),N,dt)
phi2 = compute_phi2(fired_cells_indx,N,dt,last500ms_indx,network_freq)



figure
%plot(firings(:,1),firings(:,2),'.')  %slow but works well
spy(fired_cells_indx)
title(['Raster Plot gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
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



