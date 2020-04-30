% script to make plots of gamma power vs number of cells in the PV network. 
% NOTE: these simulations do not contain PC or SOM cells.  Only PV cells. 

addpath '\\FS1\WisorData\MATLAB\Rempe\Matlab_utilities\Colormaps'  

gsyn_values = 0:0.5:6;  %0:0.05:6
Iapp_values = 200:100:900; %200:10:900
ggap_values = 0:0.5:3;
%ggap_values = 0:1:2;  % for a coarser grid
N_values    = 300:50:700;

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


% --- First figure, heat plot of gamma power as a function of the number of cells in network and gsyn
m=length(N_values);
n=length(gsyn_values);

Iapp.PV = 700;
ggap = 0;
network_freq = zeros(m,n);
phi 		 = zeros(m,n);
gamma_power1 = zeros(m,n);



parfor i=1:m
	gsyn = struct();
	
	gsyn.PVPC=0; gsyn.PVSOM=0; gsyn.PCPV=0; gsyn.PCPC=0; gsyn.PCSOM=0; gsyn.SOMPV=0; gsyn.SOMPC=0; gsyn.SOMSOM=0;
	for j=1:n
		gsyn.PVPV = gsyn_values(j);
		

		[~,avg_mem_pot,network_freq,phi,power_mat] = PV_PC_SOM_model(N_values(i),0,0,1500,0.1,gsyn,ggap,Iapp,12);
		gamma_power1(i,j) = power_mat.gamma;
		close all;
	end
end
% -----------------------------------------------------------------------------------------------------------



% --- second panel, heat plot of gamma power as a function of number of cells and Iapp ----------------------
n=length(Iapp_values);
gsyn = 3;
gamma_power2 = zeros(m,n);
parfor i=1:m
	Iapp = struct();
	Iapp.PC   = 0;
	Iapp.SOM  = 0;
	
	for j=1:n
		Iapp.PV   = Iapp_values(j);
		[~,avg_mem_pot,network_freq,phi,power_mat] = PV_PC_SOM_model(N_values(i),0,0,1500,0.1,gsyn,ggap,Iapp,12);
		gamma_power2(i,j) = power_mat.gamma;
		close all;
	end
end
% ----------------------------------------------------------------------------------------------------------


% --- third panel, heat plot of gamma power as a function of number of cells and gap conductance ----------------------
n=length(ggap_values);
gsyn = 3;
gamma_power3 = zeros(m,n);
Iapp = 700;
parfor i=1:m
	
	for j=1:n
		ggap = ggap_values(j);
		[~,avg_mem_pot,network_freq,phi,power_mat] = PV_PC_SOM_model(N_values(i),0,0,1500,0.1,gsyn,ggap,Iapp,12);
		gamma_power3(i,j) = power_mat.gamma;
		close all;
	end
end







save('gamma_heat_map_results4.23.mat','network_freq','phi','gamma_power','gsyn_values','Iapp_values')


figure
set(gcf,'Position',[100, 100, 600, 900]);
subplot(3,1,1)
imagesc(flipud(gamma_power1'))
% p1=pcolor(gamma_power1')
% shading interp
xlabel('# of cells in PV network')
ylabel('g_{syn} (nS)')
cm_plasma = plasma(100);
colormap(cm_plasma)
colorbar
title('Gamma power')
ax=gca;
ax.FontName = 'DejaVu Sans';

subplot(3,1,2)
% p2=pcolor(gamma_power2')
% shading interp
imagesc(flipud(gamma_power2'))
xlabel('# of cells in PV network')
ylabel('I_{app} (pA)')
colormap(cm_plasma)
colorbar
ax=gca;
ax.FontName = 'DejaVu Sans';

subplot(3,1,3)
% p3=pcolor(gamma_power3')
% shading interp
imagesc(flipud(gamma_power3'))
xlabel('# of cells in PV network')
ylabel('g_{gap}')
colormap(cm_plasma)
colorbar
ax=gca;
ax.FontName = 'DejaVu Sans';
