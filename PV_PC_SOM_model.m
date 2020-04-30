function [v,avg_mem_pot,network_freq,phi,firing_rate,power] = PV_PC_SOM_model(N_PV,N_PC,N_SOM,tfinal,dt,gsyn,ggap,Iapp_meanIN,Iapp_stdIN)
% 
% Usage: [v,avg_mem_potential,network_freq,phi] = PV_PC_SOM_model(N_PV,N_PC,N_SOM,tfinal,dt,gsyn,Iapp_meanIN,Iapp_stdIN)
% or     [v,avg_mem_potential,network_freq,phi] = PV_PC_SOM_model(N_PV,N_PC,N_SOM,1500,0.01) to use default values for the other inputs
% or     [v,avg_mem_potential,network_freq,phi] = PV_PC_SOM_model  to use default values for ALL parameters
% This function simulates a network made up 
% of Parvalbumin-positive (PV) cells using 
% Ferguson et al 2013, and Pyramidal cells (PC) 
% and Somatostatin cells (SOM) inspired 
% by Chen et al 2017 NEURON 1403-1418


% UPDATE 10.17.2019: default values are stored in structs gsyn, ggap,Iapp_meanIN, and Iapp_stdIN
% stored in default_param_values.mat 
% So you can call this with only 
if nargin == 0   % No arguments
  N_PV  = 500;
  N_PC  = 2000;
  N_SOM = 500;
  tfinal = 1500;
  dt = 0.01;
  load default_param_values.mat %this loads gsyn, ggap,Iapp_meanIN, and Iapp_stdIN
end 

if nargin == 5 
  load default_param_values.mat %this loads gsyn, ggap,Iapp_meanIN, and Iapp_stdIN
end 
% 
%
% INPUTS: defaults listed in brackets: 	
%     N_PV:  		    # of PV cells  [500]
%     N_PC:         # of Pyramidal cells [2000]
%     N_SOM:        # of SOM cells       [500]
% 		tfinal: 	    length of simulation (in ms) [1500]
%     dt:           timestep (in ms) [0.01]
%			gsyn: 		    synapse conductance (in nS) [1.5] (Note that this needs to be a struct with fields PVPV,PVPC,PVSOM,PCPV,PCPC,PCSOM,SOMPV,SOMPC,SOMSOM)
%			Iapp_meanIN:  mean applied current (in pA) [600]  (Also needs to be a struct)
%			Iapp_stdIN:	  std of applied current (in pA) [12] (Also needs to be a struct)
%
% OUTPUT:	v: 	last 500 ms of the voltage of all cells in the network.  
%         avg_mem_potential: membrane potential averaged across cells (contains as many elements as timesteps)

% TESTING. THIS IS USEFUL FOR RUNNING THIS FUNCTION MORE LIKE A SCRIPT
% N_PV = 500;
% N_PC = 500;
% N_SOM = 500;
% tfinal = 1500;
% dt = 0.01;
% gsyn = 1.5;
% Iapp_meanIN = 600;
% Iapp_stdIN = 12;


addpath other_models_code   % this is where compute_phi lives
addpath '\\FS1\WisorData\MATLAB\Rempe\Matlab_utilities\'

% Check if gsyn is a struct.  If not assume all values are the same and print a warning...
if ~isstruct(gsyn)
  warning('gsyn was given as a number, not a structure.  Assuming all 9 of the gsyn values are the same....')
  gsyn_in = gsyn;
  clear gsyn;
  gsyn.PVPV   = gsyn_in;  %[1.5]
  gsyn.PVPC   = gsyn_in;  %
  gsyn.PVSOM  = gsyn_in;
  gsyn.PCPV   = gsyn_in;
  gsyn.PCPC   = gsyn_in;
  gsyn.PCSOM  = gsyn_in;
  gsyn.SOMPV  = gsyn_in;
  gsyn.SOMPC  = gsyn_in;
  gsyn.SOMSOM = gsyn_in;
end

% Check to see if Iapp is a struct.  If not, assume all value are the same and print a warning..
if ~isstruct(Iapp_meanIN)
  warning('Iapp was given as a number, not a structure.  Assuming all 3 of the Iapp values are the same...')
  Iapp_mean_in = Iapp_meanIN;
  clear Iapp_meanIN;
  Iapp_meanIN.PV  = Iapp_mean_in;
  Iapp_meanIN.PC  = Iapp_mean_in;   % 0-200 seems to be a good range.  
  Iapp_meanIN.SOM = Iapp_mean_in;
end 

if ~isstruct(Iapp_stdIN)
  warning('Iapp STD was given as a number, not a structure.  Assuming all 3 of the Iapp_std values are the same... ')
  Iapp_std_in = Iapp_stdIN;
  clear Iapp_stdIN;
  Iapp_stdIN.PV  = Iapp_std_in;
  Iapp_stdIN.PC  = Iapp_std_in;
  Iapp_stdIN.SOM = Iapp_std_in;
end 



disp('Setting up.........')
tic
timesteps = tfinal/dt;  
steps_in_one_ms = 1/dt;  
t=0:dt:tfinal;
fs = (1/dt)*1000;     % sampling freq. in Hz  1000 to get from ms to sec
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
v_counts_as_fired = ones(N,1);
%esyn   = ones(N,1);

v_rest_PV  = -65;
v_rest_PC  = -60;
v_rest_SOM = -65;


% ----- PV cell parameters ---------------------------------------------------------
a(PV_cell_index)      =  0.1;         % parameter a describes time scale of recovery
b(PV_cell_index)      = -0.1;      	  % parameter b describes sensitivity of u to v
c(PV_cell_index)      = -67;          % parameter c describes after-spike reset value of v
d(PV_cell_index)      =  10;          % parameter d describes after-spike reset value of u
Cm(PV_cell_index)     =  90;
k(PV_cell_index)      =  1.7;
vr(PV_cell_index)     = -60.6;
vt(PV_cell_index)     = -43.1;        % vt=-43.1 in Ferguson et al;
k_low_PV  = 1.7;
k_high_PV = 14;
v_counts_as_fired(PV_cell_index) = 2.5; % if v reaches this value we count it as a spike
% ----------------------------------------------------------------------------------

% ----- PC cell parameters ---------------------------------------------------------
a(PC_cell_index)  = 0.02+0.08*rand(N_PC,1); %0.03;          % parameter a describes time scale of recovery
b(PC_cell_index)  = 0.2+ 0.05*rand(N_PC,1);  %-2-0.5*(rand(N_PC,1));             % parameter b describes sensitivity of u to v
c(PC_cell_index)  = -67+22*rand(N_PC,1);  %-50+15*((rand(N_PC,1)).^2);       % parameter c describes after-spike reset value of v
d(PC_cell_index)  = 8-6*((rand(N_PC,1)).^2);         % parameter d describes after-spike reset value of u was 100
Cm(PC_cell_index) = 100;
k(PC_cell_index)  = 0.7;
vr(PC_cell_index) = -60;
vt(PC_cell_index) = -40;        % vt=-43.1 in Ferguson et al;
k_low_PC  = 0.7;
k_high_PC = 0.7;
v_counts_as_fired(PC_cell_index) = 35;
% ----------------------------------------------------------------------------------

% ----- SOM cell parameters ---------------------------------------------------------
a(SOM_cell_index)  =  0.001;      % parameter a describes time scale of recovery
b(SOM_cell_index)  = 0.2;         % parameter b describes sensitivity of u to v
c(SOM_cell_index)  = -64;         % parameter c describes after-spike reset value of v
d(SOM_cell_index)  =  8;          % parameter d describes after-spike reset value of u
Cm(SOM_cell_index) =  100;
k(SOM_cell_index)  =  1.7;
vr(SOM_cell_index) = -64;
vt(SOM_cell_index) = -35;         % vt=-43.1 in Ferguson et al;
k_low_SOM  = 0.34;
k_high_SOM = 1.5;
v_counts_as_fired(SOM_cell_index) = 11;  
% ------------------------------------------------------------------------------------

% --- Connectivity parameters --------------------------------------------------------
PVPV_conn_prob   = 0.12;
PVPC_conn_prob   = 0.9;  % probably 1 according to Pfeffer et al 2013?  
PVSOM_conn_prob  = 0;
PCPV_conn_prob   = 0;
PCPC_conn_prob   = 0.5;
PCSOM_conn_prob  = 0;
SOMPV_conn_prob  = 0;
SOMPC_conn_prob  = 0;
SOMSOM_conn_prob = 0;  % should be 0, according to Chen et al NEURON 2017 and 
                          % Gibson et al Nature 1999 (calls them LTS neurons)
% ------------------------------------------------------------------------------------



% ----- Synaptic parameters (may want to make these specific to cell type too) -------
alpha = (1/0.27);
beta  = (1/1.8);
s_inf = alpha/(alpha+beta);
tau_s = 1/(alpha+beta);
esyn  = -85;   %*ones(N,1);

% -------------------------------------------------------------------------------------

% ----- Gap junction parameters -------------------------------------------------------
if ~isstruct(ggap)
  g_gap_SOM = ggap; %1 %2.1 +- 2.3 nS according to Gibson 1999 recorded in rat ;  % in nS   
  g_gap_PV  = ggap;
else
  g_gap_SOM = ggap.SOM;
  g_gap_PV  = ggap.PV;
end 
% -------------------------------------------------------------------------------------

% ----  applied current (all other synaptic inputs, following Ferguson et al) ---------
% Iapp_std  = Iapp_stdIN*ones(N,1);
% Iapp_mean = Iapp_meanIN*ones(N,1);
% Iapp  = randn(N,1).*Iapp_std + Iapp_mean;  % Set this up once, don't update it, like Scott said (Frances' postdoc)
Iapp_std  = zeros(N,1);
Iapp_mean = zeros(N,1);
Iapp_mean(PV_cell_index)  = Iapp_meanIN.PV;
Iapp_mean(PC_cell_index)  = Iapp_meanIN.PC;
Iapp_mean(SOM_cell_index) = Iapp_meanIN.SOM;
Iapp_std(PV_cell_index)   = Iapp_stdIN.PV;
Iapp_std(PC_cell_index)   = Iapp_stdIN.PC;
Iapp_std(SOM_cell_index)  = Iapp_stdIN.SOM; 
Iapp  = randn(N,1).*Iapp_std + Iapp_mean;
% -------------------------------------------------------------------------------------


%  ------------ SYNAPTIC CONNECTIVITY ------------------------------------------------------------
conn_mat = build_connection_matrix(N_PV,N_PC,N_SOM,PV_cell_index,PC_cell_index,SOM_cell_index,...
                      PVPV_conn_prob,PVPC_conn_prob,PVSOM_conn_prob, ...
                      PCPV_conn_prob,PCPC_conn_prob,PCSOM_conn_prob, ...
                      SOMPV_conn_prob,SOMPC_conn_prob,SOMSOM_conn_prob, ...
                      gsyn.PVPV,gsyn.PVPC,gsyn.PVSOM,gsyn.PCPV,gsyn.PCPC,gsyn.PCSOM,gsyn.SOMPV,gsyn.SOMPC,gsyn.SOMSOM);
disp(unique(conn_mat))
mapping = cell(N,1);          % cell array listing all presynaptic cells for each post-syn cell
for i=1:N
  mapping{i} = find(conn_mat(i,:));  %this just has locations.  Make it include values too somehow,
end
conductance = cell(N,1);
for row=1:N
  for i=1:length(mapping{row})
    conductance{row}(i) = conn_mat(row,mapping{row}(i));
  end
end                               % conductance contains the values of all the gsyn, not just locations.  
cells_receiving_syn_input = find(~cellfun(@isempty,mapping));
% -----------------------------------------------------------------------------------------------

% ------------ GAP JUNCTION CONNECTIVITY FOR SOM POPULATION -------------------------------------
gap_conn_mat_SOM = rand(N_SOM);
locs_low     = find(gap_conn_mat_SOM<=0.93);  % 93% connectivity between SOM cells
locs_high    = find(gap_conn_mat_SOM> 0.93);
gap_conn_mat_SOM(locs_low)  = 1;
gap_conn_mat_SOM(locs_high) = 0;
gap_conn_mat_SOM = gap_conn_mat_SOM-triu(gap_conn_mat_SOM)+(tril(gap_conn_mat_SOM))';  % make it symmetric
gap_conn_mat_SOM = gap_conn_mat_SOM.*(1-eye(N_SOM));                                % remove main diagonal
%gap_conn_mat = g_gap*gap_conn_mat;
% -----------------------------------------------------------------------------------------------

% ------------ GAP JUNCTION CONNECTIVITY FOR PV POPULATION -------------------------------------
gap_conn_mat_PV = rand(N_PV);
locs_low     = find(gap_conn_mat_PV<=0.02);  % 2% connectivity between PV cells? (see Bartos 2002)
locs_high    = find(gap_conn_mat_PV> 0.02);
gap_conn_mat_PV(locs_low)  = 1;
gap_conn_mat_PV(locs_high) = 0;
gap_conn_mat_PV = gap_conn_mat_PV-triu(gap_conn_mat_PV)+(tril(gap_conn_mat_PV))';  % make it symmetric
gap_conn_mat_PV = gap_conn_mat_PV.*(1-eye(N_PV));                                % remove main diagonal
%gap_conn_mat = g_gap*gap_conn_mat;
% -----------------------------------------------------------------------------------------------


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
Isyn = zeros(N,1);			% initialize vector for synaptic input
Igap = zeros(N,1);      % initialize vector for gap junction input
% ---------------------------------------------------------------------------------------




disp('Done initializing and setting up.  ')
disp('Now updating the network....')
tic
% Bootstrap: Take one step of FE before starting Adams-Bashforth method (works as of 3.28.2019)
% s_now = s(:,1);
% CS = zeros(N,1);
% for i=1:N
%   CS(i) = sum(s_now(mapping{i}));
% end

s_now = s(:,1);
CS = zeros(N,1);
for i=1:length(cells_receiving_syn_input) 
  cond_this_cell = conductance{cells_receiving_syn_input(i)};   % this should give a row vector of gsyn values for all cells presynaptic to cell i
  CS(i) = cond_this_cell*(s_now(mapping{cells_receiving_syn_input(i)}));  % dot product: row vector multiplied by a column vec 
end


%Isyn = (gsyn*CS).*(v(:,1)-esyn); % Old way, when gsyn was constant for whole network
Isyn = (CS).*(v(:,1)-esyn);

uold = u;
Isyn_old = Isyn;
k_old = k;
v(:,2) = v(:,1)+(dt./Cm).*(k.*(v(:,1)-vr).*(v(:,1)-vt)-u+Iapp-Isyn);
unew = u+dt.*a.*(b.*(v(:,1)-vr)-u); 
vplot(1)           = v(1,1);
v_with_spikes(:,1) = v(:,1);    
uplot(1)           = u(1);







% ---------------------------- UPDATING LOOP ---------------------------------------------
% ----------------------------------------------------------------------------------------
for step=2:timesteps           
  vplot(step)           = v(1,step);
  v_with_spikes(:,step) = v(:,step);        % Before you reset v, save it
  uplot(step)           = u(1);
  
	fired = find(v(:,step)>=v_counts_as_fired);    			 % indices of cells that spiked
  fired_cells_indx(:,step) = v(:,step)>=v_counts_as_fired;            % logical array of whether cell fire
  %firings  = [firings; t(step)+0*fired,fired]; % also works, but is slow
  

   v_high_PV    = find(v(PV_cell_index,step)>vt(PV_cell_index));    % for the PV cells, find which ones are above threshold
   v_low_PV     = find(v(PV_cell_index,step)<=vt(PV_cell_index));   % and set k for these cells. 
   k(v_high_PV) = k_high_PV;                                        
   k(v_low_PV)  = k_low_PV;
   v_high_PC    = find(v(PC_cell_index,step)>vt(PC_cell_index));    % for the PC cells, find which ones are above threshold
   v_low_PC     = find(v(PC_cell_index,step)<=vt(PC_cell_index));   % and set k for these cells.
   k(v_high_PC+N_PV) = k_high_PC;                                        
   k(v_low_PC+N_PV)  = k_low_PC;
   v_high_SOM   = find(v(SOM_cell_index,step)>vt(SOM_cell_index));  % for the SOM cells, find which ones are above threshold
   v_low_SOM    = find(v(SOM_cell_index,step)<=vt(SOM_cell_index)); % and set k for these cells .
   k(v_high_SOM+N_PV+N_PC) = k_high_SOM;                                       
   k(v_low_SOM+N_PV+N_PC)  = k_low_SOM;

  % v_high_all = find(v(:,step)>vt(:));  % for the PV cells, find which ones are above threshold
  % v_low_all  = find(v(:,step)<=vt(:));
  % k(v_high_all) = k_high_PV;                          % set k just for PV cells, the other types may not need this
  % k(v_low_all)  = k_low_PV;  
    %k(v(fired,step)>=vt(fired)) = k_high(fired);  %Faster. logical indexing
    %k(v(fired,step)<vt(fired))  = k_low(fired); 





  if ~isempty(fired)                        % if at least one cell fired this timestep, reset v,u,T, and k for those cells
    v(fired,step) = c(fired);               % after spike reset v
    u(fired)      = u(fired)+d(fired);	          % after spike reset u
    T(fired,step:step+steps_in_one_ms) = 1;
  end
  
  % --- synapses ----------------------------
  % Old way, works great, but is a bit slow
  %Isyn(:,step) = (gsyn*conn_mat*s(:,step)).*(v(:,step)-esyn);
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
  if g_gap_SOM ~= 0
    vdiffs_SOM = gap_conn_mat_SOM.*(repmat(v(SOM_cell_index,step)',N_SOM,1) - repmat(v(SOM_cell_index,step),1,N_SOM));
    Igap(SOM_cell_index) = g_gap_SOM*sum(vdiffs_SOM,2);
  end 
  % -----------------------------------------------------

  % --- gap junctions for PV cells ----------------------
  if g_gap_PV ~= 0
    vdiffs_PV = gap_conn_mat_PV.*(repmat(v(PV_cell_index,step)',N_PV,1) - repmat(v(PV_cell_index,step),1,N_PV));
    Igap(PV_cell_index) = g_gap_PV*sum(vdiffs_PV,2);
  end
  % -----------------------------------------------------


  % Forward Euler (works 3.22.2019)
  % v(:,step+1) = v(:,step)+(dt./Cm).*(k.*(v(:,step)-vr).*(v(:,step)-vt)-u+Iapp-Isyn(:,step)-Igap); % this is Forward Euler.  Do something better. 
  % u = u+dt.*a.*(b.*(v(:,step)-vr)-u); 
  
  % Adams-Bashforth (works as of 3.22.2019)
  % v(:,step+1) = v(:,step) + (dt/2).*(3*fv(v(:,step),Cm,k,vr,vt,u,Isyn(:,step),Iapp)-fv(v(:,step-1),Cm,k_old,vr,vt,uold,Isyn(:,step-1),Iapp));
  % unew        = u         + (dt/2).*(3*fu(u,v(:,step),a,b,vr)-fu(uold,v(:,step-1),a,b,vr));
  % uold = u;
  % u = unew;
  % %Isyn_old = Isyn;
  % k_old = k;

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

  % Old way, too slow
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
    disp(['25% finished.... Should be finished at ',num2str(datestr(now + 3*seconds(et))), ' which is ',num2str(4*et), ' seconds.'])
  elseif step == round(timesteps/2)
    disp('50% finished...')
  elseif step == round(0.75*timesteps)
    disp('75% finished...')
  end 

end  % end of looping through timesteps
disp('Done updating the model. Starting post-processing...')


avg_mem_pot.all = mean(v_with_spikes,1);
avg_mem_pot.PV  = mean(v_with_spikes(PV_cell_index,:),1);
avg_mem_pot.PC  = mean(v_with_spikes(PC_cell_index,:),1);
avg_mem_pot.SOM = mean(v_with_spikes(SOM_cell_index,:),1);


% --- filter the signals here before computing frequency info -----------------------------------
% hpFilt = designfilt('highpassfir','StopbandFrequency',0.001,...
%                     'PassbandFrequency',0.01,'PassbandRipple',100,...
%                     'StopbandAttenuation',25,'DesignMethod','kaiserwin');  % This one was from JW
hpFilt = designfilt('highpassiir', ...
         'StopBandFrequency',0.1,'PassbandFrequency',0.5,'StopbandAttenuation',5,'PassbandRipple',0.2, ...
         'SampleRate',fs,'DesignMethod','ellip','MatchExactly','stopband');

% avg_mem_pot.all = filtfilt(hpFilt,avg_mem_pot.all);
% avg_mem_pot.PV  = filtfilt(hpFilt,avg_mem_pot.PV);
% avg_mem_pot.PC  = filtfilt(hpFilt,avg_mem_pot.PC);
% avg_mem_pot.SOM = filtfilt(hpFilt,avg_mem_pot.SOM);
% -----------------------------------------------------------------------------------------------



% ----- compute network frequency following Ferguson et al 2013 --------------------------------- 
% (freq at which there is a spectral peak in the last 500ms of population activity)
signal_all = avg_mem_pot.all(last_500ms_indices)-mean(avg_mem_pot.all(last_500ms_indices));  % remove mean to detrend
if N_PV > 0 signal_PV  = avg_mem_pot.PV(last_500ms_indices) -mean(avg_mem_pot.PV(last_500ms_indices)); ,end  % remove mean to detrend 
if N_PC > 0 signal_PC  = avg_mem_pot.PC(last_500ms_indices) -mean(avg_mem_pot.PC(last_500ms_indices)); ,end   % remove mean to detrend
if N_SOM> 0 signal_SOM = avg_mem_pot.SOM(last_500ms_indices)-mean(avg_mem_pot.SOM(last_500ms_indices));,end  % remove mean to detrend
% -----------------------------------------------------------------------------------------------



[pxx_all,f_all] = pwelch(signal_all,fs/2,round(0.6*fs/2),2*fs,fs);
if N_PV > 0 [pxx_PV,f_PV]   = pwelch(signal_PV,fs/2,round(0.6*fs/2),2*fs,fs); ,end 
if N_PC > 0 [pxx_PC,f_PC]   = pwelch(signal_PC,fs/2,round(0.6*fs/2),2*fs,fs); ,end 
if N_SOM> 0 [pxx_SOM,f_SOM] = pwelch(signal_SOM,fs/2,round(0.6*fs/2),2*fs,fs);,end 

cutoff_f_all = find(f_all==250);  % To find the max frequency, only consider frequencies below 250 Hz
if N_PV>0  cutoff_f_PV  = find(f_PV ==250);,end
if N_PC>0  cutoff_f_PC  = find(f_PC ==250);,end 
if N_SOM>0 cutoff_f_SOM = find(f_SOM==250);,end 

max_freq_indx_all = find(pxx_all==max(pxx_all(1:cutoff_f_all)));
if N_PV>0  max_freq_indx_PV  = find(pxx_PV ==max(pxx_PV(1:cutoff_f_PV)));   ,end 
if N_PC>0  max_freq_indx_PC  = find(pxx_PC ==max(pxx_PC(1:cutoff_f_PC)));   ,end 
if N_SOM>0 max_freq_indx_SOM = find(pxx_SOM==max(pxx_SOM(1:cutoff_f_SOM))); ,end 

network_freq.all = f_all(max_freq_indx_all);
if N_PV>0  network_freq.PV  = f_PV(max_freq_indx_PV);  ,end
if N_PC>0  network_freq.PC  = f_PC(max_freq_indx_PC);  ,end 
if N_SOM>0 network_freq.SOM = f_SOM(max_freq_indx_SOM);,end 
% ------------------------------------------------------------------------------------------------


% ----------- Compute power in delta, theta, gamma, and alpha bands ------------------------------
if N_PV>0  power.PV  = compute_power_in_freq_bands(avg_mem_pot.PV,dt,last_500ms_indices); ,end
if N_PC>0  power.PC  = compute_power_in_freq_bands(avg_mem_pot.PC,dt,last_500ms_indices); ,end 
if N_SOM>0 power.SOM = compute_power_in_freq_bands(avg_mem_pot.SOM,dt,last_500ms_indices);,end 

if N_PV>0  power.PV.gamma  ,end 
if N_PC>0  power.PC.gamma  ,end 
if N_SOM>0 power.SOM.gamma ,end 
% TimeFreqMatPV  = mf_tfcm(signal_PV, 6,[1:200],fs,1,0.21,'power');   % Use Morlet wavelet to view power in freq ranges vs time
% TimeFreqMatPC  = mf_tfcm(signal_PC, 6,[1:200],fs,1,0.21,'power');   % Use Morlet wavelet to view power in freq ranges vs time
% TimeFreqMatSOM = mf_tfcm(signal_SOM,6,[1:200],fs,1,0.21,'power');  % Use Morlet wavelet to view power in freq ranges vs time
if N_PV>0  [modulusPV,phasesPV,TransformPV] = morlet_wav(signal_PV,fs,1/(5*pi),1,200,1);    ,end %1/(5*pi) is a nice balance between freq and temp resolution
if N_PC>0  [modulusPC,phasesPC,TransformPC] = morlet_wav(signal_PC,fs,1/(5*pi),1,200,1);    ,end  %1/(5*pi) is a nice balance between freq and temp resolution
if N_SOM>0 [modulusSOM,phasesSOM,TransformSOM] = morlet_wav(signal_SOM,fs,1/(5*pi),1,200,1);,end   %1/(5*pi) is a nice balance between freq and temp resolution


% -- These lines only work if all 3 cell types are present ----------------------------------------
% [cfsPV,f_cwtPV]   = cwt(signal_PV,'amor',fs);    % matlab code for Morlet wavelwt using the wavelet toolbox
% [cfsPC,f_cwtPC]   = cwt(signal_PV,'amor',fs);    % matlab code for Morlet wavelwt using the wavelet toolbox
% [cfsSOM,f_cwtSOM] = cwt(signal_PV,'amor',fs);    % matlab code for Morlet wavelwt using the wavelet toolbox

% cwt_argsPV  = {t(last_500ms_indices),f_cwtPV,abs(cfsPV).^2};
% cwt_argsPC  = {t(last_500ms_indices),f_cwtPC,abs(cfsPC).^2};
% cwt_argsSOM = {t(last_500ms_indices),f_cwtSOM,abs(cfsSOM).^2};
% ------------------------------------------------------------------------------------------------

% ---- Compute average firing rates in last 0.5 seconds ------------------------------------------
save('test_fired_cells_indx.mat','fired_cells_indx')
if N_PV>0  firing_rate.PV = 2*mean(sum(fired_cells_indx(PV_cell_index,last_500ms_indices),2));  ,end % the 2* is because I'm using 0.5 seconds of data
if N_PC>0  firing_rate.PC = 2*mean(sum(fired_cells_indx(PC_cell_index,last_500ms_indices),2));  ,end % the 2* is because I'm using 0.5 seconds of data
if N_SOM>0 firing_rate.SOM = 2*mean(sum(fired_cells_indx(SOM_cell_index,last_500ms_indices),2));,end  % the 2* is because I'm using 0.5 seconds of data
% ------------------------------------------------------------------------------------------------


% ------------ Compute network coherency ---------------------------------------------------------
if N_PV>0  phi.PV  = compute_phi2(fired_cells_indx(PV_cell_index,:),N_PV,dt,last500ms_indx,network_freq.PV);   ,end 
if N_PC>0  phi.PC  = compute_phi2(fired_cells_indx(PC_cell_index,:),N_PC,dt,last500ms_indx,network_freq.PC);   ,end 
if N_SOM>0 phi.SOM = compute_phi2(fired_cells_indx(SOM_cell_index,:),N_SOM,dt,last500ms_indx,network_freq.SOM);,end 
% ------------------------------------------------------------------------------------------------


% --------- PLOTTING ------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
toc
gsyn
figure
%plot(firings(:,1),firings(:,2),'.')  % slow but works when you uncomment firings line above
spy(fired_cells_indx)
title(['Raster Plot gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
ylabel('Cell number')
xlabel('Timesteps')
axis normal
axis([timesteps-50/dt timesteps 0 N])

figure
plot(t(1:end-1),v_with_spikes(1,:))
title(['Output of cell 1, gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Voltage (mV)')


figure
plot(t(1:end-1),splot)
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Synaptic gating variable s'  )

figure
plot(t(1:end-1),avg_mem_pot.all)
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, All cells)')

if N_PV>0
figure
plot(t(1:end-1),avg_mem_pot.PV)
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, PV)')
end 

if N_PC>0
figure
plot(t(1:end-1),avg_mem_pot.PC)
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, PC)')
end

if N_SOM>0
figure
plot(t(1:end-1),avg_mem_pot.SOM)
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Time (ms)')
ylabel('Population activity (average membrane potential, SOM)')
end 

figure
plot(f_all,10*log10(pxx_all),f_all,10*log10(pxx_all),'x')
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) All types of cells')
axis([0 350 ylim])

if N_PV>0
figure
plot(f_PV,10*log10(pxx_PV),f_PV,10*log10(pxx_PV),'x')
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) PV cells')
axis([0 350 ylim])
end 

if N_PC>0
figure
plot(f_PC,10*log10(pxx_PC),f_PC,10*log10(pxx_PC),'x')
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) PC cells')
axis([0 350 ylim])
end 

if N_SOM>0
figure
plot(f_SOM,10*log10(pxx_SOM),f_SOM,10*log10(pxx_SOM),'x')
title(['gsyn.PVPV = ', num2str(gsyn.PVPV), ' Iapp = ', num2str(Iapp_mean(1))])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB/Hz) SOM cells')
axis([0 350 ylim])
end 


if N_PV>0
figure
pcolor(modulusPV)
shading interp
colorbar
xlabel('Time ')
ylabel('Frequency (Hz)')
title('PV cells only')
end 

if N_PC>0
figure
pcolor(modulusPC)
shading interp
colorbar
xlabel('Time ')
ylim([1 200])
ylabel('Frequency (Hz)')
title('PC cells only')
end 

if N_SOM>0
figure
pcolor(modulusSOM)
shading interp
colorbar
xlabel('Time ')
ylim([1 200])
ylabel('Frequency (Hz)')
title('SOM cells only')
end 

% figure
% surf(cwt_argsPV{:},'edgecolor','none');
% view(0,90);
% axis tight;
% shading interp; colormap(parula(128));
% h=colorbar;
% h.Label.String = 'Power';
% xlabel('Time')
% ylabel('Hz')


