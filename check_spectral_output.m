% script to check spectral output of rodent sleep data
% This script works and produces expected output as of July 13, 2018.  
% read in edf file and accompanying edf file
addpath 'C:\Users\wisorlab\Documents\MATLAB\mrempe\Epoch-Based-Processing\Timed_Intervals\export\AutomatedScoring\AutoScore_with_GUI\AutoScore_using_edf\New_version_multiple_sleep_states_edf'

inputFile = '\\FS1\WisorData\Rempe_Data\Data\test_data\E2697Base2.edf'; %E2697Base.edf 

% uncomment next two lines if you need to load the data.  commented to save time once data were loaded
[fileHeader,channelHeader,data,state,fs,textdata,epoch_length] = LoadAndMergeEdfAndTxt_MJR(inputFile);
data = data{1};







% --- using periodogram       ---
% --- first with all the data ---
freq   =  fs(1);
window = hamming(length(data));
[pxxAll,fxAll] = periodogram(data,window,length(data),freq,'power');
figure
plot(fxAll,pxxAll)
title('Experimental data all states using periodogram')
ylabel('')

% SWS only 
figure
S_indices = find(strcmp(state(:),'S'));
window = hamming(length(S_indices));
[pxxS,fxS] = periodogram(data(S_indices),window,length(S_indices),freq,'power');
plotS=plot(fxS,pxxS,'r');
ax=gca;
ax.XLim = [0 100];
%title('Costa output SWS using periodogram')
hold on 

% REMS only 
R_indices = find(strcmp(state(:),'R'));
window = hamming(length(R_indices));
[pxxR,fxR] = periodogram(data(R_indices),window,length(R_indices),freq,'power');
plotR=plot(fxR,pxxR,'Color',[1 0.5 0]);
ax=gca;
ax.XLim = [0 100];
%title('Costa output REMS using periodogram')
plotR.Color(4) = 0.2;

% Wake only 
W_indices = find(strcmp(state(:),'W'));
window = hamming(length(W_indices));
[pxxW,fxW] = periodogram(data(W_indices),window,length(W_indices),freq,'power');
plotW=plot(fxW,pxxW);
ax=gca;
ax.XLim = [0 100];
title('Experimental data using periodogram')
hold off
legend('SWS','REMS','Wake')
plotW.Color(4) = 0.2;  

% --- now compare to using pwelch ---
% pwelch was much slower and seemed to give less helpful results so I'm commenting it out.  7.5.2018
% --- first with all the data -------
% disp('Now trying pwelch....')
% [pxx,f] = pwelch(data,200,[],1024,freq);
% figure
% plot(f,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('Power Spectral Density (dB/Hz)')
% title('Experimental data all states using pwelch')


% % pwelch on SWS
% [pxxS,fS]=pwelch(data(S_indices),200,[],1024,freq);
% figure
% plot(fS,10*log10(pxxS))
% xlabel('Frequency (Hz)')
% ylabel('Power Spectral Density (dB/Hz)')
% title('Experimental data SWS using pwelch')


% % pwelch on wake 
% [pxxW,fW]=pwelch(data(W_indices),200,[],1024,freq);
% figure
% plot(fW,10*log10(pxxW))
% xlabel('Frequency (Hz)')
% ylabel('Power Spectral Density (dB/Hz)')
% title('Experimental data Wake using pwelch')

% % pwelch on REMS 
% [pxxR,fR]=pwelch(data(R_indices),200,[],1024,freq);
% figure
% plot(fR,10*log10(pxxR))
% xlabel('Frequency (Hz)')
% ylabel('Power Spectral Density (dB/Hz)')
% title('Experimental data REMS using pwelch')









