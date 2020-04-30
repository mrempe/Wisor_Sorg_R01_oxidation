% script to run PV_model3.m many times to make a figure like Figure 3 in Ferguson et al 2013,
% except this time making a heat map of gamma power at each gsyn/Iapp combination, rather than phi

addpath '\\FS1\WisorData\MATLAB\Rempe\Matlab_utilities\Colormaps'  % where the new colormaps live

gsyn = 0:0.2:8;
Iapp = 200:20:1000;


network_freq = zeros(length(gsyn),length(Iapp));
phi = zeros(length(gsyn),length(Iapp));


m=length(gsyn);
n=length(Iapp);

parfor i=1:m
	for j=1:n
		[~,~,network_freq(i,j),phi(i,j),spec_pwr] = PV_model3(500,1500,0.01,gsyn(i),Iapp(j),12);
		gamma(i,j) = spec_pwr.gamma;
		theta(i,j) = spec_pwr.theta;
		close all;
	end
end

save('gI_heat_map_results3.7.19.mat','network_freq','phi','gamma','theta')





figure
imagesc(flipud(gamma'))
xlabel('gsyn')
ylabel('Iapp')
cm_plasma = plasma(100);
colormap(cm_plasma)
title('Power in gamma')
% uncomment the next lines to put the correct tick labels on the axes.  
ax = gca;
% ax.XTick = 0:length(gsyn);
% ax.XTickLabel = {'0', '0', };  % Add the values you need, making sure XTickLabel has the same number of elements as XTick
% ax.YTick = 1:2:length(Iapp);
% ax.YTickLabel = {'900', '850', '800', '750', '700', ...};  % Add values you need, making sure to go backwards so small values are at bottom