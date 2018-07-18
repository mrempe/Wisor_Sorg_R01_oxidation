function dydt = Costa(t,y)
% .m file to specify the differential equations in the Cost et al 2016 model
%

% --- All the parameter values ---
% --------------------------------
% ---                          ---
% --- Sleep Model parameters   ---
% Membrane time in [s]  
tau_W 	= 1500E3;
tau_N 	= 600E3;
tau_R 	= 60E3;

% Neurotransmitter time constants in [ms] 
tau_E 	= 25E3;
tau_G 	= 10E3;
tau_A 	= 10E3;

% Maximum firing rate in [s^-1] */
F_W_max	= 6.5;
F_N_max	= 5.;
F_R_max	= 5.;

% Sigmoid slope parameters in [aU] */
alpha_W	= 0.5;
alpha_N	= 0.175;
alpha_R	= 0.13;

% Sigmoid threshold parameters in [aU] */
beta_W		= -0.4;
beta_R		= -0.9;

% Neurotransmitter release scaling in [s^-1] */
gamma_E		= 5.;
gamma_G		= 4.;
gamma_A		= 2.;

% Weights for neurotransmitter efficacy in [aU]*/
g_GW		= -1.68;
g_AW		= 1.;
g_GR		= -1.3;
g_AR		= 1.6;
g_ER		= -4.;
g_EN		= -2.;

% Sleep Homeostasis parameter */
H_max		= 1.;		% in [aU] */
theta_W		= 2.;		% in [s] */
tau_hw		= 34830E3;	% 580.5 min in [s] */
tau_hs		= 30600E3;	% 510 min in [s] */
kappa		= 1.5;		% in [aU] */



% set up one large vector y, and dydt

% ordering of variables from Cost paper to one large vector y
%
% me 		Costa
% y(1)		F_W			% sleep regulatory network
% y(2) 		F_N
% y(3) 		F_R
% y(4) 		C_E
% y(5) 		C_G
% y(6)		C_A
% y(7)		h
%
% y(8)		V_p
% y(9)		V_i
% y(10)		s_ep
% y(11)		s_gp
% y(12)		s_ei
% y(13)		s_gi
% y(14)		[Na]
% y(15)		g_KNa
% y(16)		sigma_p


dydt = zeros(16,1);

dydt(1) = (F_W_max *0.5*(1 + tanh((g_GW * y(5) + g_AW * y(6) - beta_W)/alpha_W)))/tau_W;
dydt(2) = (F_N_max *0.5*(1 + tanh((g_EN * y(4) + kappa*y(7))/alpha_N)))/tau_N;
dydt(3) = (F_R_max *0.5*(1 + tanh((g_ER * y(4) + g_GR * C_G[N] + g_AR * C_A[N] - beta_R)/alpha_R))     - f_R [N])/tau_R;
dydt(4) = (tanh(y(1)/gamma_E) - y(4))/tau_E;
dydt(5) = (tanh(y(2)/gamma_G) - y(5))/tau_G;
dydt(6) = (tanh(y(3)/gamma_A) - y(6))/tau_A;
dydt(7) = (H_max-y(7))/tau_hw*Heaviside(f_W[N]-theta_W) - h[N]/tau_hs*Heaviside(theta_W- f_W[N]);

