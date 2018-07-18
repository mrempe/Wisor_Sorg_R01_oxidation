function make_shaded_state_percentages_plot(t,sleep_state,sleep_dep_start_stop_times)

% this needs to read in several runs of the model and plot the average percentages in each
% 2 hour window

sleep_dep_start_times = sleep_dep_start_stop_times(1:2:end-1);
sleep_dep_end_times	  = sleep_dep_start_stop_times(2:2:end);


window_length_in_hours = 2;  % Janne uses 2 hour windows in the figure in the paper

window_length_in_elements = find(abs(t-window_length_in_hours)<1e-15);
window_times 			  = t(ceil(window_length_in_elements/2)):window_length_in_hours:t(end-ceil(window_length_in_elements/2));


for j=1:size(sleep_state,2)
	index = 1;
	for i=1:window_length_in_elements:length(t)-window_length_in_elements
		window_indices 			 = i:i+window_length_in_elements -1;
		wake_percentage(index,j) = length(find(sleep_state(window_indices,j)=='W'))/window_length_in_elements;
		SWS_percentage(index,j)  = length(find(sleep_state(window_indices,j)=='S'))/window_length_in_elements;
		REM_percentage(index,j)  = length(find(sleep_state(window_indices,j)=='R'))/window_length_in_elements;
		index = index+1;
	end
end

% Now average over the trials
wake_percentage_means = mean(wake_percentage,2)';
SWS_percentage_means  = mean(SWS_percentage,2)';
REM_percentage_means  = mean(REM_percentage,2)';




FigHandle = figure;
set(FigHandle,'Position',[300 50 941, 216]);

plot(window_times,REM_percentage_means,window_times,REM_percentage_means+SWS_percentage_means)
yl=ylim;
hold on 
Y = [REM_percentage_means' SWS_percentage_means' wake_percentage_means' ];
h = area(window_times,Y,'LineStyle','none');
h(1).FaceColor = [0.2863 0 0.5725];       % purple
h(2).FaceColor = [0.4275 0.7137 1]; %light blue
h(3).FaceColor = [0.8588 0.4275 0];  %orange



if ~isempty(sleep_dep_start_stop_times)
	gray = [0.31 0.31 0.31];
	for i=1:length(sleep_dep_start_times)
		p = patch([sleep_dep_start_times(i) sleep_dep_end_times(i) sleep_dep_end_times(i) sleep_dep_start_times(i)],[0 0 1 1],[0.5725 0.2863 0]);
		p.EdgeColor='none';
		
	end 
end 

hold off 
set(gca,'XTick',0:24:t(end));	% set major tick marks at 24 hour intervals
set(gca,'YTick',0:0.1:1);
if sleep_dep_start_times(1)==26
	V=axis;
	axis([0 120 V(3) V(4)])
end 
ha=gca;
ha.XAxis.MinorTick = 'on';
ha.XAxis.MinorTickValues = ha.XAxis.Limits(1):8:ha.XAxis.Limits(2);  % set minor tick marks at 8 hour intervals
ha.TickDir    = 'out';
ha.TickLength = [0.02 0.02];
ha.LineWidth  = 1.5;
ha.FontSize   = 12;
ha.YTick 	  = 0:0.2:1;
ha.YTickLabel = {'0%','20%','40%','60%','80%','100%'};
set(gca,'box','off')
set(gca,'color','none')
ylabel('% of total time')
ylim([0 1.1])
% Add the yellow and black bars at the top to indicate light/dark cycle
%rectangle('Position',[xlowerleft ylowerleft width height ])
hold on
rectangle('Position',[0 1 12 0.1],'EdgeColor','none','FaceColor','y')
rectangle('Position',[12 1 12 0.1],'EdgeColor','none','FaceColor','k')
rectangle('Position',[24 1 12 0.1],'EdgeColor','none','FaceColor','y')
rectangle('Position',[36 1 12 0.1],'EdgeColor','none','FaceColor','k')
rectangle('Position',[48 1 12 0.1],'EdgeColor','none','FaceColor','y')
rectangle('Position',[60 1 12 0.1],'EdgeColor','none','FaceColor','k')
rectangle('Position',[72 1 12 0.1],'EdgeColor','none','FaceColor','y')
rectangle('Position',[84 1 12 0.1],'EdgeColor','none','FaceColor','k')
rectangle('Position',[96 1 12 0.1],'EdgeColor','none','FaceColor','y')
rectangle('Position',[108 1 12 0.1],'EdgeColor','none','FaceColor','k')
if t(end)>120
	rectangle('Position',[120 1 12 0.1],'EdgeColor','none','FaceColor','y')
end
hold off

% Add the labels of Baseline 'Day1', 'Day2', etc.  
text(8.5,-0.08,{'Baseline'},'FontSize',12)
text(34,-0.08,{'Day1'},'FontSize',12)
text(58,-0.08,{'Day2'},'FontSize',12)
text(82,-0.08,{'Day3'},'FontSize',12)
text(106,-0.08,{'Day4'},'FontSize',12)
if t(end)>120
	text(127,-0.08,{'Day5'},'FontSize',12)
end
