% script to plot gamma power and freq with most power using only PC cells.  

addpath '\\FS1\WisorData\MATLAB\Rempe\Matlab_utilities\Colormaps'  

gsyn_values = 0:0.5:6;  %0:0.05:6
Iapp_values = 0:100:900; %200:10:900
ggap_values = 0;
%ggap_values = 0:1:2;  % for a coarser grid
N_values    = 1000:100:2200;

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
	
	gsyn.PVPC=0; gsyn.PVSOM=0; gsyn.PCPV=0; gsyn.PVPV=0; gsyn.PCSOM=0; gsyn.SOMPV=0; gsyn.SOMPC=0; gsyn.SOMSOM=0;
	for j=1:n
		gsyn.PCPC = gsyn_values(j);
		

		[~,avg_mem_pot,network_freq,phi,power_mat] = PV_PC_SOM_model(0,N_values(i),0,1500,0.1,gsyn,ggap,Iapp,12);
		gamma_power1(i,j) = power_mat.PC.gamma;
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
	Iapp.PV   = 0;
	Iapp.SOM  = 0;
	
	for j=1:n
		Iapp.PC   = Iapp_values(j);
		[~,avg_mem_pot,network_freq,phi,power_mat] = PV_PC_SOM_model(0,N_values(i),0,1500,0.1,gsyn,ggap,Iapp,12);
		gamma_power2(i,j) = power_mat.PC.gamma;
		close all;
	end
end
% ----------------------------------------------------------------------------------------------------------


% Finally, do a plot of peak frequency as a function of Iapp
peak_freq = zeros(n,1);
for j=1:n
	Iapp.PC   = Iapp_values(j);
	[~,avg_mem_pot,network_freq,phi,power_mat] = PV_PC_SOM_model(0,2000,0,1500,0.1,gsyn,ggap,Iapp,12);
	peak_freq(j) = network_freq.PC;
	close all;
end










save('PC_only_gamma_heat_map_results5.2.mat','network_freq','phi','gamma_power1','gamma_power2','peak_freq','gsyn_values','Iapp_values')


figure
set(gcf,'Position',[100, 100, 600, 900]);
subplot(3,1,1)
imagesc(flipud(gamma_power1'))
% p1=pcolor(gamma_power1')
% shading interp
xlabel('# of cells in PC network')
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
xlabel('# of cells in PC network')
ylabel('I_{app} (pA)')
colormap(cm_plasma)
colorbar
ax=gca;
ax.FontName = 'DejaVu Sans';

subplot(3,1,3)
p=plot(Iapp_values,peak_freq,'x-')
p.LineWidth = 2
ax=gca;
ax.FontName = 'DejaVu Sans';
set(gca,'Color','none')
xlabel('Iapp')
ylabel('Peak Frequency')

