function F = F_Costa_onGPU(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19)
% .m file to specify the differential equations in the Costa et al 2016 model
% for use with the Chang_Math_Comp_1987 algorithm
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
F_W_max	= 6.5e-3;
F_N_max	= 5.0e-3;
F_R_max	= 5.0e-3;

% Sigmoid slope parameters in [aU] */
alpha_W	= 0.5;
alpha_N	= 0.175;
alpha_R	= 0.13;

% Sigmoid threshold parameters in [aU] */
beta_W	= -0.4;
beta_R	= -0.9;

% Neurotransmitter release scaling in [s^-1] */
gamma_E	= 5.0e-3;
gamma_G	= 4.0e-3;
gamma_A	= 2.0e-3;

% Weights for neurotransmitter efficacy in [aU]*/
g_GW	= -1.68;
g_AW	= 1.0;
g_GR	= -1.3;
g_AR	= 1.6;
g_ER	= -4.0;
g_EN	= -2.0;

% Sleep Homeostasis parameter */
H_max	= 1.;		% in [aU] */
theta_W	= 2e-3;		% in [ms^-1] */
tau_hw	= 34830E3;	% 580.5 min in [s] */
tau_hs	= 30600E3;	% 510 min in [s] */
kappa	= 1.5;		% in [aU] */

% --- Cortex Model parameters   ---
% Membrane time in ms */
tau_p 		= 30;
tau_i 		= 30;

% Maximum firing rate in ms^-1 */
Qp_max		= 30e-3;
Qi_max		= 60e-3;

% Sigmoid threshold in mV */
theta_p		= -58.5;
theta_i		= -58.5;

% Sigmoid gain in mV */
sigma_p_0	= 7;
tau_s		= 100;
sigma_i		= 6;

% Scaling parameter for sigmoidal mapping (dimensionless) */
C1			= pi/sqrt(3);

% Parameters of the firing adaption */
alpha_Na	= 2.0;			% Sodium influx per spike  in mM ms 	*/
tau_Na		= 1.0;			% Sodium time constant	    in ms		*/

R_pump   	= 0.09;        	% Na-K pump constant	    in mM/ms 	*/
Na_eq    	= 9.5;         	% Na-eq concentration	    in mM		*/

% PSP rise time in ms^-1 */
gamma_e		= 70E-3;
gamma_g		= 58.6E-3;

% Membrane capacitance */
C_m			= 1;

% Leak weight in aU */
g_L    		= 1.;

% Synaptic weight in ms */
g_AMPA 		= 1.0;
g_GABA 		= 1.0;

% KNa condutivity in mS/m^2 */
g_KNa_0		= 1.33;
tau_g		= 10;

% Reversal potentials in mV */
% Synaptic */
E_AMPA  	= 0;
E_GABA  	= -70;

% Leak */
E_L_p 		= -66;
E_L_i 		= -64;

% Potassium */
E_K    		= -100;

% Noise parameters in ms^-1 */
mphi		= 0.0;
dphi		= 2;
input		= 0.0;

% Connectivities (dimensionless) */
% Label indicates N_{from -> to} */
N_pp		= 120;
N_ip		= 72;
N_pi		= 90;
N_ii		= 90;
% -------------------------------------
% --- END OF PARAMETERS ---------------



% set up one large vector y, and dydt

% ordering of variables from Costa paper to one large vector y
%
% me 		Costa
% F(1)		F_W			% sleep regulatory network
% F(2) 		F_N
% F(3) 		F_R
% F(4) 		C_E
% F(5) 		C_G
% F(6)		C_A
% F(7)		h
%
% F(8)		V_p			% Cortex model
% F(9)		V_i
% F(10)		s_ep
% F(11)		s_ei
% F(12)		s_gp
% F(13)		s_gi
% F(14)		x_ep (first derivative of s_ep)
% F(15) 	x_ei (first derivative of s_ei)
% F(16) 	x_gp (first derivative of s_gp)
% F(17) 	x_gi (first derivative of x_gi)
% F(18)		[Na]
% F(19)		g_KNa
% F(20)		sigma_p
	



F(1) = (F_W_max *0.5*(1 + tanh((g_GW * X5 + g_AW * X6 - beta_W)/alpha_W)) - X1)/tau_W;
F(2) = (F_N_max *0.5*(1 + tanh((g_EN * X4 + kappa*X7)/alpha_N))           - X2)/tau_N;
F(3) = (F_R_max *0.5*(1 + tanh((g_ER * X4 + g_GR *X5+ g_AR * X6 - beta_R)/alpha_R))-X3)/tau_R;
F(4) = (tanh(X1)/gamma_E) - X4)/tau_E;
F(5) = (tanh(X2)/gamma_G) - X5)/tau_G;
F(6) = (tanh(X3)/gamma_A) - X6)/tau_A;
F(7) = ((H_max-X7)/tau_hw)*heaviside(X1-theta_W) - (X7/tau_hs)*heaviside(theta_W-X1);

% cortex model
w_KNa = 0.37/(1+(38.7/X18)^3.5);
I_KNa = X19 * w_KNa * (X8 - E_K);
F(8)  = -(g_L * (X8-E_L_p) + g_AMPA * X1) * (X8 - E_AMPA) + g_GABA * X12 * (X8 - E_GABA) )/tau_p - (1/C_m) * I_KNa;
F(9)  = -(g_L * (X9-E_L_i) + g_AMPA * X1) * (X9 - E_AMPA) + g_GABA * X13 * (X9 - E_GABA) )/tau_i;
F(10) = X14;
F(11) = X15;
F(12) = X16;
F(13) = X17;
get_Qp = Qp_max / (1 + exp(-C1 * (X8 - theta_p) / X20));
get_Qi = Qi_max / (1 + exp(-C1 * (X9 - theta_i) / sigma_i));
Na_pump = R_pump*(X18^3/(X18^3+3375) - Na_eq^3/(Na_eq^3+3375));
% dydt(14) = gamma_e^2 * (N_pp * get_Qp +dphi*(0.1/(4*sqrt(3)))*randn - X(10)) - 2 * gamma_e * X(14);  %  dphi multiplies randn because randn gives SD of 1
% dydt(15) = gamma_e^2 * (N_ip * get_Qp +dphi*(0.1/(4*sqrt(3)))*randn - X(11)) - 2 * gamma_e * X(15);
F(14) = gamma_e^2 * (N_pp * get_Qp - X10) - 2 * gamma_e * X14;  %  dphi multiplies randn because randn gives SD of 1
F(15) = gamma_e^2 * (N_ip * get_Qp - X11) - 2 * gamma_e * X15;
F(16) = gamma_g^2 * (N_pi * get_Qi - X12) - 2 * gamma_g * X16;
F(17) = gamma_g^2 * (N_ii * get_Qi - X13) - 2 * gamma_g * X17;
F(18) = (alpha_Na * get_Qp - Na_pump)/tau_Na;
F(19) = (g_KNa_0 * (2*X5)*(1-0.6*X4)*(1-0.95*X6) - X19)/tau_g;
F(20) = (sigma_p_0 - (4*X4 + 2*X6) - X20)/tau_s;


F=F';  % make it a column vector, not a row vector