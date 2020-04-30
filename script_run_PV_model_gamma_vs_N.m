% Script to run PV_model3.m with different numbers of cells in the network to see how sensitive gamma is to the # of cells in the
% network

N = 100:10:800;


network_freq = zeros(length(gsyn),length(Iapp));
phi = zeros(length(gsyn),length(Iapp));


m=length(N);
gamma = zeros(size(N));
theta = zeros(size(N));



parfor i=1:m
	[~,~,network_freq(i),phi(i),spec_pwr] = PV_model3(N(i),1500,0.01,3.8,580,12);
	gamma(i) = spec_pwr.gamma;
	theta(i) = spec_pwr.theta;
	close all;
end

save('gamma_vs_N_results3.8.19.mat','network_freq','phi','gamma','theta')

figure
p=plot(N,gamma,N,gamma,'.');
p(2).MarkerSize = 16;

