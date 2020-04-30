% script to run PV_PC_SOM_model.m many times to make a figure like Figure 3 in Ferguson et al 2013, except showing 
% power in the gamma frequency range rather than coherence

addpath '\\FS1\WisorData\MATLAB\Rempe\Matlab_utilities\Colormaps'  % where the new colormaps live

tstart = tic;

gsyn_values = 0:0.5:6;  %0:0.05:6 		% commented values were used on 4.23
Iapp_values = 200:100:1000; %200:10:900	% commented values were used on 4.23

Iapp=struct();
Iapp.PV  = 0;
Iapp.PC  = 0;
Iapp.SOM = 0;

gsyn = struct();
gsyn.PVPV = 0;
gsyn.PVPC = 0;
gsyn.PVSOM = 0;
gsyn.PCPV = 0;
gsyn.PCPC = 0;
gsyn.PCSOM = 0;
gsyn.SOMPV = 0;
gsyn.SOMPC = 0;
gsyn.SOMSOM = 0;

m=length(gsyn_values);
n=length(Iapp_values);


network_freq = zeros(m,n);
phi 		 = zeros(m,n);
gamma_power  = zeros(m,n);

%gamma_power = parallel.pool.Constant(gamma_power);
gsyn_values = parallel.pool.Constant(gsyn_values);
Iapp_values = parallel.pool.Constant(Iapp_values);

parfor i=1:m
	gsyn = struct();
	Iapp = struct();
	Iapp.PV   = 0;
	Iapp.SOM  = 0;
	gsyn.PVPC=0; gsyn.PVSOM=0; gsyn.PCPV=0; gsyn.PVPV=0; gsyn.PCSOM=0; gsyn.SOMPV=0; gsyn.SOMPC=0; gsyn.SOMSOM=0;
	for j=1:n
		gsyn.PCPC = gsyn_values.Value(i);
		Iapp.PC   = Iapp_values.Value(j);

		[~,avg_mem_pot,network_freq,phi_mat,power_mat] = PV_PC_SOM_model(0,500,0,1500,0.05,gsyn,0,Iapp,12);
		gamma_power(i,j) = power_mat.PC.gamma;
		phi(i,j)         = phi_mat.PC;
		close all;
	end
end

%save('gamma_heat_map_results4.22.mat','network_freq','phi','gamma_power','gsyn_values','Iapp_values')

toc(tstart)

figure
%imagesc(flipud(gamma_power'))
pcolor(gamma_power')
shading interp
xlabel('g_{syn} (nS)')
ylabel('I_{app} (pA)')
cm_plasma = plasma(100);
colormap(cm_plasma)
title('Gamma power')
% % uncomment the next lines to put the correct tick labels on the axes.  
% ax = gca;

figure
pcolor(phi')
shading interp
xlabel('g_{syn} (nS)')
ylabel('I_{app} (pA)')
cm_plasma = plasma(100);
colormap(cm_plasma)
title('phi')
