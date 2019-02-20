% script to run PV_model2.m many times to make a figure like Figure 3 in Ferguson et al 2013

gsyn = 0:0.25:6;
Iapp = 400:25:900;


network_freq = zeros(length(gsyn),length(Iapp));
phi = zeros(length(gsyn),length(Iapp));


m=length(gsyn);
n=length(Iapp);

parfor i=1:m
	for j=1:n
		[~,~,network_freq(i,j),phi(i,j)] = PV_model2(500,1500,0.01,gsyn(i),Iapp(j),12);
		close all;
	end
end

save('Fig3_repro_results2.19.mat','network_freq','phi')


figure
imagesc(flipud(phi'))
xlabel('gsyn')
ylabel('Iapp')
colormap('gray')
% uncomment the next lines to put the correct tick labels on the axes.  
ax = gca;
% ax.XTick = 0:length(gsyn);
% ax.XTickLabel = {'0', '0', };  % Add the values you need, making sure XTickLabel has the same number of elements as XTick
% ax.YTick = 1:2:length(Iapp);
% ax.YTickLabel = {'900', '850', '800', '750', '700', ...};  % Add values you need, making sure to go backwards so small values are at bottom