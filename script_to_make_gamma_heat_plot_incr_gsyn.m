% script to make a heat plot of how increasing GABA (gsyn) changes gamma power in the PV network
% NOTE:  This is for a network of 500 PV cells.  No other types of cells.  


addpath '\\FS1\WisorData\MATLAB\Rempe\Matlab_utilities\Colormaps'



% load the gamma power data (as function of gsyn and Iapp)
load 'BEST_gamma_heat_map_gsyn_Iapp_fine_res4.22.mat'


% or re-run the code to generate that data
%script_runmodel_gamma_heatmap



% ----- MAKE Figure ----------------------------------------
% ----------------------------------------------------------
% --- Panel 1: regular heat plot of gamma power ------------
figure
set(gcf,'Position',[100, 100, 600, 900]);
s1=subplot(2,1,1)
pcolor(gamma_power')
shading interp
box off
xlabel('g_{syn} (nS)')
ylabel('I_{app} (pA)')
cm_plasma = plasma(100);
colormap(cm_plasma)
colorbar
title('Power in gamma range (80-100 Hz) ')
ax = gca;
ax.YTick=[1 10 20 30 40 50 60 70];
ax.YTickLabel={'200', '300', '400', '500','600','700', '800','900'};
ax.XTick = [1 20:20:120];
ax.XTickLabel={'0' '1' '2' '3' '4' '5' '6'};
ax.FontName = 'DejaVu Sans';
s1.FontSize = 12;
% -----------------------------------------------------------

% smallest possible is gsyn_values(2)-gsyn_values(1)
amount_to_increase_gsyn = 1;  % nS
num_panels_in_increase = fix(amount_to_increase_gsyn/(gsyn_values(2)-gsyn_values(1)));
N=length(gsyn_values);

% ---- Panel 2: heat plot of the change in gamma that results from changing GABA
gamma_flipped = gamma_power';
gamma_diffs   = zeros(size(gamma_flipped,1),N-num_panels_in_increase);
for j=1:N-num_panels_in_increase
	gamma_diffs(:,j) = gamma_flipped(:,j+num_panels_in_increase)-gamma_flipped(:,j);
end

s2=subplot(2,1,2)
pcolor(gamma_diffs)
shading interp
box off
xlabel('g_{syn} (nS)')
ylabel('I_{app} (pA)')
colormap(cm_plasma)
colorbar
title('Change in gamma power due to increasing GABA')
ax=gca;
ax.YTick=[1 10 20 30 40 50 60 70];
ax.YTickLabel={'200', '300', '400', '500','600','700', '800','900'};
ax.XTick = [1 20:20:120];
ax.XTickLabel={'0' '1' '2' '3' '4' '5' '6'};
ax.FontName = 'DejaVu Sans';
s2.FontSize = 12;
