function [v,avg_mem_pot_all,avg_mem_pot_PV,avg_mem_pot_PC,avg_mem_pot_SOM,network_freq_all,network_freq_PV,network_freq_PC,network_freq_SOM,phi] = PV_PC_SOM_model(N_PV,N_PC,N_SOM,tfinal,dt,gsyn,Iapp_meanIN,Iapp_stdIN)
% 
% Usage: [v,avg_mem_potential,network_freq,phi] = PV_PC_SOM_model(N_PV,N_PC,N_SOM,tfinal,dt,gsyn,Iapp_meanIN,Iapp_stdIN)
%
% This function simulates a network made up 
% of Parvalbumin-positive (PV) cells using 
% Ferguson et al 2013, and Pyramidal cells (PC) 
% and Somatostatin cells (SOM) inspired 
% by Chen et al 2017 NEURON 1403-1418

%
% INPUTS: defaults listed in brackets: 	
%     N_PV:  		    # of PV cells  [500]
%     N_PC:         # of Pyramidal cells [?]
%     N_SOM:        # of SOM cells       [?]
% 		tfinal: 	    length of simulation (in ms) [1500]
%     dt:           timestep (in ms) [0.01]
%			gsyn: 		    synapse conductance (in nS) [150 NOTE typo in Ferguson et al 2013]
%			Iapp_meanIN:  mean applied current (in pA) [600]
%			Iapp_stdIN:	  std of applied current (in pA) [12]
%
% OUTPUT:	v: 	last 500 ms of the voltage of all cells in the network.  
%         avg_mem_potential: membrane potential averaged across cells (contains as many elements as timesteps)

% TESTING. THIS IS USEFUL FOR RUNNING THIS FUNCTION MORE LIKE A SCRIPT
% N_PV = 500;
% N_PC = 200;
% N_SOM = 200;
% tfinal = 1500;
% dt = 0.01;
% gsyn = 150;
% Iapp_meanIN = 600;
% Iapp_stdIN = 12;
profile on

addpath other_models_code   % this is where compute_phi lives
addpath '\\FS1\WisorData\MATLAB\Rempe\Matlab_utilities'  % this is where progressbar.m lives


disp('Setting up.........')
tic
timesteps = tfinal/dt;  
steps_in_one_ms = 1/dt;  
t=0:dt:tfinal;
last500ms_indx = find(abs(t-(tfinal-500))<1e-12);  %ONLY ONE Element: the location of t vector where t=500ms from end
last_500ms_indices = find(t>=tfinal-500);
last_500ms_indices = last_500ms_indices(1:end-1);  % remove the last one so sizes work

% index the cells based on type
N = N_PV + N_PC + N_SOM;           % total # of cells
PV_cell_index  = 1:N_PV;           % the first N_PV cells are PV
PC_cell_index  = N_PV+1:N_PV+N_PC; % next N_PC cells are PC
SOM_cell_index = N_PV+N_PC+1:N;    % next N_SOM cells are SOM

% set up all parameter vectors first
a      = ones(N,1);        
b      = ones(N,1);        
c      = ones(N,1);        
d      = ones(N,1);        
S      = rand(N,N);  
Cm     = ones(N,1);
k      = ones(N,1);
vr     = ones(N,1);
vt     = ones(N,1);
k_low  = ones(N,1);
k_high = ones(N,1);

v_rest_PV  = -65;
v_rest_PC  = -65;
v_rest_SOM = -65;


% ----- PV cell parameters ---------------------------------------------------------
a(PV_cell_index)      =  0.1;         % parameter a describes time scale of recovery
b(PV_cell_index)      = -0.1;      	  % parameter b describes sensitivity of u to v
c(PV_cell_index)      = -70;          % parameter c describes after-spike reset value of v
d(PV_cell_index)      =  10;          % parameter d describes after-spike reset value of u
S(PV_cell_index)      =  0.5;  
Cm(PV_cell_index)     =  90;
k(PV_cell_index)      =  1.7;
vr(PV_cell_index)     = -60.6;
vt(PV_cell_index)     = -43.1;        % vt=-43.1 in Ferguson et al;
k_low_PV  = 1.7;
k_high_PV = 14;
% ----------------------------------------------------------------------------------

% ----- PC cell parameters ---------------------------------------------------------
a(PC_cell_index)  =  0.1;         % parameter a describes time scale of recovery
b(PC_cell_index)  = -0.1;         % parameter b describes sensitivity of u to v
c(PC_cell_index)  = -70;          % parameter c describes after-spike reset value of v
d(PC_cell_index)  =  10;          % parameter d describes after-spike reset value of u
S(PC_cell_index)  =  0.5;  
Cm(PC_cell_index) =  90;
k(PC_cell_index)  =  1.7;
vr(PC_cell_index) = -60.6;
vt(PC_cell_index) = -43.1;        % vt=-43.1 in Ferguson et al;
k_low_PC  = 1.7;
k_high_PC = 14;
% ----------------------------------------------------------------------------------

% ----- SOM cell parameters ---------------------------------------------------------
a(SOM_cell_index)  =  0.1;         % parameter a describes time scale of recovery
b(SOM_cell_index)  = -0.1;         % parameter b describes sensitivity of u to v
c(SOM_cell_index)  = -70;          % parameter c describes after-spike reset value of v
d(SOM_cell_index)  =  10;          % parameter d describes after-spike reset value of u
S(SOM_cell_index)  =  0.5;  
Cm(SOM_cell_index) =  90;
k(SOM_cell_index)  =  1.7;
vr(SOM_cell_index) = -60.6;
vt(SOM_cell_index) = -43.1;        % vt=-43.1 in Ferguson et al;
k_low_SOM  = 1.7;
k_high_SOM = 14;
% ------------------------------------------------------------------------------------

% --- Connectivity parameters --------------------------------------------------------
PVPV_conn_prob   = 0.12;
PVPC_conn_prob   = 0.12;
PVSOM_conn_prob  = 0.12;
PCPV_conn_prob   = 0.12;
PCPC_conn_prob   = 0.12;
PCSOM_conn_prob  = 0.12;
SOMPV_conn_prob  = 0.12;
SOMPC_conn_prob  = 0.12;
SOMSOM_conn_prob = 0;
% ------------------------------------------------------------------------------------



% ----- Synaptic parameters (may want to make these specific to cell type too) -------
alpha = (1/0.27);
beta  = (1/1.8);
s_inf = alpha/(alpha+beta);
tau_s = 1/(alpha+beta);
esyn  = -85;   %*ones(N,1);   
% -------------------------------------------------------------------------------------


% ----  applied current (all other synaptic inputs, following Ferguson et al) ---------
Iapp_std  = Iapp_stdIN*ones(N,1);
Iapp_mean = Iapp_meanIN*ones(N,1);
Iapp  = randn(N,1).*Iapp_std + Iapp_mean;  % Set this up once, don't update it, like Scott said (Frances' postdoc)
% -------------------------------------------------------------------------------------


%  ------------ CONNECTIVITY ------------------------------------------------------------
conn_mat = build_connection_matrix(N_PV,N_PC,N_SOM,PV_cell_index,PC_cell_index,SOM_cell_index,...
                      PVPV_conn_prob,PVPC_conn_prob,PVSOM_conn_prob, ...
                      PCPV_conn_prob,PCPC_conn_prob,PCSOM_conn_prob, ...
                      SOMPV_conn_prob,SOMPC_conn_prob,SOMSOM_conn_prob);


% [post_syn_cells,pre_syn_cells] = find(conn_mat);
% post_syn_cells = unique(post_syn_cells);
% pre_syn_cells  = unique(pre_syn_cells);

% for i=post_syn_cells'
%   pre_syn_connections{i} = find(conn_mat(i,:));  % cell array where each entry corresponds to a post-synaptic cell. Each entry lists pre-synaptic cells, if any, that connect to this cell
% end

% % ---- PV connectivity matrix, connect each PV cell randomly --------------------------- 
% % ---- with all other PV cells with probability 12% ------------------------------------
% PV_conn_mat  = rand(N_PV);    % PV connectivity matrix for PV-PV connections.  Columns are pre-synaptic cells, rows are post-synaptic
% PV_locs_low  = find(PV_conn_mat <= 0.12);
% PV_locs_high = find(PV_conn_mat > 0.12);
% PV_conn_mat(PV_locs_low)  = 1;
% PV_conn_mat(PV_locs_high) = 0;
% [PV_post_syn_cells,PV_pre_syn_cells] = find(PV_conn_mat);
% PV_post_syn_cells = unique(PV_post_syn_cells);
% PV_pre_syn_cells  = unique(PV_pre_syn_cells);
% for i=PV_post_syn_cells'
% 	PV_pre_syn_connections{i} = find(PV_conn_mat(i,:));  % cell array where each entry corresponds to a post-synaptic cell. Each entry lists pre-synaptic cells, if any, that connect to this cell
% end
% % ---------------------------------------------------------------------------------------

% % ---- PC connectivity matrix, connect each PC cell randomly --------------------------- 
% % ---- with all other PC cells with probability 12% ------------------------------------
% PC_conn_mat  = rand(N_PC);    % PV connectivity matrix for PC-PC connections.  Columns are pre-synaptic cells, rows are post-synaptic
% PC_locs_low  = find(PC_conn_mat <= 0.12);
% PC_locs_high = find(PC_conn_mat > 0.12);
% PC_conn_mat(PC_locs_low)  = 1;
% PC_conn_mat(PC_locs_high) = 0;
% [PC_post_syn_cells,PC_pre_syn_cells] = find(PC_conn_mat);
% PC_post_syn_cells = unique(PC_post_syn_cells);
% PC_pre_syn_cells  = unique(PC_pre_syn_cells);
% for i=PC_post_syn_cells'
%   PC_pre_syn_connections{i} = find(PC_conn_mat(i,:));  % cell array where each entry corresponds to a post-synaptic cell. Each entry lists pre-synaptic cells, if any, that connect to this cell
% end

% % ----------------------------------------------------------------------------------------


% % ---- SOM connectivity matrix, connect each SOM cell randomly --------------------------- 
% % ---- with all other SOM cells with probability 12% -------------------------------------
% SOM_conn_mat  = rand(N_SOM);    % PV connectivity matrix for SOM-SOM connections.  Columns are pre-synaptic cells, rows are post-synaptic
% SOM_locs_low  = find(SOM_conn_mat <= 0.12);
% SOM_locs_high = find(SOM_conn_mat > 0.12);
% SOM_conn_mat(SOM_locs_low)  = 1;
% SOM_conn_mat(SOM_locs_high) = 0;
% [SOM_post_syn_cells,SOM_pre_syn_cells] = find(SOM_conn_mat);
% SOM_post_syn_cells = unique(SOM_post_syn_cells);
% SOM_pre_syn_cells  = unique(SOM_pre_syn_cells);
% for i=SOM_post_syn_cells'
%   SOM_pre_syn_connections{i} = find(SOM_conn_mat(i,:));  % cell array where each entry corresponds to a post-synaptic cell. Each entry lists pre-synaptic cells, if any, that connect to this cell
% end
% % ---------------------------------------------------------------------------------------

% % ---- concatenate the lists of presynaptic cells and post-synaptic, being careful to offset because of different cell types
% shifted_pre_connections_PC   = cellfun(@(x) x+N_PV,PC_pre_syn_connections,'un',0);
% shifted_pre_connections_SOM  = cellfun(@(x) x+N_PV+N_PC,SOM_pre_syn_connections,'un',0);
% %shifted_post_connections_PC  = cellfun(@(x) x+N_PV,PC_pre_syn_connections,'un',0);
% %shifted_post_connections_SOM = cellfun(@(x) x+N_PV+N_PC,SOM_pre_syn_connections,'un',0); 
% pre_syn_connections = [PV_pre_syn_connections shifted_pre_connections_PC  shifted_pre_connections_SOM];
% post_syn_cells      = [PV_post_syn_cells;      PC_post_syn_cells+N_PV; SOM_post_syn_cells+N_PV+N_PC];



disp('Initializing...............')
% ----- Initialize variables (v,u,s) ----------------------------------------------------
v_PV  = 10*rand(N_PV,1) + v_rest_PV*ones(N_PV,1); %-65*ones(N,1);    % Initial values of v, all -65 mV
v_PC  = 10*rand(N_PC,1) + v_rest_PC*ones(N_PC,1);
v_SOM = 10*rand(N_SOM,1)+ v_rest_SOM*ones(N_SOM,1);

v=[v_PV;v_PC;v_SOM];
u = b.*v;            % Initial values of u done in original Izhevich code
%u = zeros(size(v));  % try initializing u to zero
s = zeros(length(v),timesteps);
firings=[];             % spike timings
T    = zeros(N,timesteps);		  % boolean matrix that contains a 1 in row i if 
                                % cell i has fired in the last ms. rows correspond to cell numbers, columns to time steps
Isyn = zeros(N,timesteps);			% initialize vector for synaptic input
% ---------------------------------------------------------------------------------------




disp('Done initializing and setting up.  ')
disp('Now updating the network....')
% ---------------------------- UPDATING LOOP ---------------------------------------------
% ----------------------------------------------------------------------------------------
for step=1:timesteps           
  vplot(step)           = v(1,step);
  v_with_spikes(:,step) = v(:,step);        % Before you reset v, save it
  uplot(step)           = u(1);
  
	fired = find(v(:,step)>=2.5);    			 % indices of cells that spiked
  fired_cells_indx(:,step) = v(:,step)>=2.5;            % logical array of whether cell fire
  %firings  = [firings; t(step)+0*fired,fired]; % also works, but is slow
  

  v_high_PV = find(v(PV_cell_index,step)>vt(PV_cell_index));  % for the PV cells, find which ones are above threshold
  v_low_PV  = find(v(PV_cell_index,step)<=vt(PV_cell_index));
  k(v_high_PV) = k_high_PV;                          % set k just for PV cells, the other types may not need this
  k(v_low_PV)  = k_low_PV;
    %k(v(fired,step)>=vt(fired)) = k_high(fired);  %Faster. logical indexing
    %k(v(fired,step)<vt(fired))  = k_low(fired); 





  if ~isempty(fired)                        % if at least one cell fired this timestep, reset v,u,T, and k for those cells
    v(fired,step) = c(fired);               % after spike reset v
    u(fired)      = u(fired)+d(fired);	          % after spike reset u
    T(fired,step:step+steps_in_one_ms) = 1;
  end
  
  Isyn(:,step) = (gsyn*conn_mat*s(:,step)).*(v(:,step)-esyn);


  v(:,step+1) = v(:,step)+(dt./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn(:,step)); % this is Forward Euler.  Do something better. 
  u = u+dt.*a.*(b.*(v(:,step)-vr)-u); 
  
  for i=1:N % loop over cells
    if T(i,step)
      s(i,step+1) = s_inf + (s(i,step)-s_inf)*exp(-dt/tau_s);
    else
      s(i,step+1) = s(i,step)*exp(-beta*dt);
    end
  end



  splot(step) = s(5,step);
  
  if step ==round(timesteps/4)
    disp('25% finished...')
  elseif step == round(timesteps/2)
    disp('50% finished...')
  elseif step == round(0.75*timesteps)
    disp('25% finished')
  end 

end  % end of looping through timesteps
disp('Done updating the model. Starting post-processing')


avg_mem_pot_all = mean(v_with_spikes,1);
avg_mem_pot_PV  = mean(v_with_spikes(PV_cell_index,:),1);
avg_mem_pot_PC  = mean(v_with_spikes(PC_cell_index,:),1);
avg_mem_pot_SOM = mean(v_with_spikes(SOM_cell_index,:),1);


% ----- compute network frequency following Ferguson et al 2013 --------------------------------- 
% CHANGE THIS TO BE CELL SPECIFIC: 
% (freq at which there is a spectral peak in the last 500ms of population activity)
fs = (1/dt)*1000;     % sampling freq. in Hz  1000 to get from ms to sec
signal_all = avg_mem_pot_all(100001:end)-mean(avg_mem_pot_all(100001:end));  % remove mean to detrend
signal_PV  = avg_mem_pot_PV(100001:end) -mean(avg_mem_pot_PV(100001:end));  % remove mean to detrend
signal_PC  = avg_mem_pot_PC(100001:end) -mean(avg_mem_pot_PC(100001:end));  % remove mean to detrend
signal_SOM = avg_mem_pot_SOM(100001:end)-mean(avg_mem_pot_SOM(100001:end));  % remove mean to detrend

[pxx_all,f_all] = pwelch(signal_all,fs/2,round(0.6*fs/2),200000,fs);
[pxx_PV,f_PV]   = pwelch(signal_PV,fs/2,round(0.6*fs/2),200000,fs);
[pxx_PC,f_PC]   = pwelch(signal_PC,fs/2,round(0.6*fs/2),200000,fs);
[pxx_SOM,f_SOM] = pwelch(signal_SOM,fs/2,round(0.6*fs/2),200000,fs);

max_freq_indx_all = find(pxx_all==max(pxx_all));
max_freq_indx_PV  = find(pxx_PV==max(pxx_PV));
max_freq_indx_PC  = find(pxx_PC==max(pxx_PC));
max_freq_indx_SOM = find(pxx_SOM==max(pxx_SOM));

network_freq_all  = f_all(max_freq_indx_all)
network_freq_PV   = f_PV(max_freq_indx_PV)
network_freq_PC   = f_PC(max_freq_indx_PC)
network_freq_SOM  = f_SOM(max_freq_indx_SOM)
% ------------------------------------------------------------------------------------------------


% ------------ Compute network coherency ---------------------------------------------------------
phi.PV  = compute_phi2(fired_cells_indx(PV_cell_index,last_500ms_indices),N_PV,dt,last500ms_indx,network_freq_PV);
phi.PC  = compute_phi2(fired_cells_indx(PC_cell_index,last_500ms_indices),N_PC,dt,last500ms_indx,network_freq_PC);
phi.SOM = compute_phi2(fired_cells_indx(SOM_cell_index,last_500ms_indices),N_SOM,dt,last500ms_indx,network_freq_SOM);
% ------------------------------------------------------------------------------------------------


% --------- PLOTTING ------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
toc
figure
%plot(firings(:,1),firings(:,2),'.')  % slow but works when you uncomment firings line above
spy(fired_cells_indx)
title(['Raster Plot gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(250))])
ylabel('Cell number')
xlabel('Timesteps')
axis normal
axis([timesteps-50/dt timesteps 0 N])

figure
plot(t(1:end-1),v_with_spikes(1,:))
title(['Output of cell 1, gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Voltage (mV)')


figure
plot(t(1:end-1),splot)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Synaptic gating variable s'  )

figure
plot(t(1:end-1),avg_mem_pot_all)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, All cells)')


figure
plot(t(1:end-1),avg_mem_pot_PV)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, PV)')

figure
plot(t(1:end-1),avg_mem_pot_PC)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, PC)')

figure
plot(t(1:end-1),avg_mem_pot_SOM)
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, SOM)')


figure
plot(f_all,10*log10(pxx_all),f_all,10*log10(pxx_all),'x')
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) All types of cells')
axis([0 350 ylim])

figure
plot(f_PV,10*log10(pxx_PV),f_PV,10*log10(pxx_PV),'x')
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) PV cells')
axis([0 350 ylim])

figure
plot(f_PC,10*log10(pxx_PC),f_PC,10*log10(pxx_PC),'x')
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) PC cells')
axis([0 350 ylim])

figure
plot(f_SOM,10*log10(pxx_SOM),f_SOM,10*log10(pxx_SOM),'x')
title(['gsyn = ', num2str(gsyn), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) SOM cells')
axis([0 350 ylim])

profile viewer

