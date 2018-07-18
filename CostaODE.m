function dydt = CostaODE(t,y)
% .m file to specify the differential equations in the Costa et al 2016 model
%

% --- All the parameter values ---
% --------------------------------
% ---                          ---
% --- Sleep Model parameters   ---
% Membrane time in [ms]  
tau_W 	= 1500E3;
tau_N 	= 600E3;
tau_R 	= 60E3;

% Neurotransmitter time constants in [ms] 
tau_E 	= 25E3;
tau_G 	= 10E3;
tau_A 	= 10E3;

% Maximum firing rate in [ms^-1] */
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

% Neurotransmitter release scaling in [ms^-1] */
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
% y(1)		F_W			% sleep regulatory network
% y(2) 		F_N
% y(3) 		F_R
% y(4) 		C_E
% y(5) 		C_G
% y(6)		C_A
% y(7)		h
%
% y(8)		V_p			% Cortex model
% y(9)		V_i
% y(10)		s_ep
% y(11)		s_ei
% y(12)		s_gp
% y(13)		s_gi
% y(14)		x_ep (first derivative of s_ep)
% y(15) 	x_ei (first derivative of s_ei)
% y(16) 	x_gp (first derivative of s_gp)
% y(17) 	x_gi (first derivative of x_gi)
% y(18)		[Na]
% y(19)		g_KNa
% y(20)		sigma_p
	


dydt = zeros(20,1);

dydt(1) = (F_W_max *0.5*(1 + tanh((g_GW * y(5) + g_AW * y(6) - beta_W)/alpha_W)) - y(1))/tau_W;
dydt(2) = (F_N_max *0.5*(1 + tanh((g_EN * y(4) + kappa*y(7))/alpha_N))           - y(2))/tau_N;
dydt(3) = (F_R_max *0.5*(1 + tanh((g_ER * y(4) + g_GR *y(5)+ g_AR * y(6) - beta_R)/alpha_R))-y(3))/tau_R;
dydt(4) = (tanh(y(1)/gamma_E) - y(4))/tau_E;
dydt(5) = (tanh(y(2)/gamma_G) - y(5))/tau_G;
dydt(6) = (tanh(y(3)/gamma_A) - y(6))/tau_A;
dydt(7) = ((H_max-y(7))/tau_hw)*heaviside(y(1)-theta_W) - (y(7)/tau_hs)*heaviside(theta_W-y(1));

% cortex model
w_KNa = 0.37/(1+(38.7/y(18))^3.5);
I_KNa = y(19) * w_KNa * (y(8) - E_K);
dydt(8)  = -(g_L * (y(8)-E_L_p) + g_AMPA * y(10) * (y(8) - E_AMPA) + g_GABA * y(12) * (y(8) - E_GABA) )/tau_p - (1/C_m) * I_KNa;
dydt(9)  = -(g_L * (y(9)-E_L_i) + g_AMPA * y(11) * (y(9) - E_AMPA) + g_GABA * y(13) * (y(9) - E_GABA) )/tau_i;
dydt(10) = y(14);
dydt(11) = y(15);
dydt(12) = y(16);
dydt(13) = y(17);
get_Qp = Qp_max / (1 + exp(-C1 * (y(8) - theta_p) / y(20)));
get_Qi = Qi_max / (1 + exp(-C1 * (y(9) - theta_i) / sigma_i));
Na_pump = R_pump*(y(18)^3/(y(18)^3+3375) - Na_eq^3/(Na_eq^3+3375));
% dydt(14) = gamma_e^2 * (N_pp * get_Qp +dphi*(0.1/(4*sqrt(3)))*randn - y(10)) - 2 * gamma_e * y(14);  %  dphi multiplies randn because randn gives SD of 1
% dydt(15) = gamma_e^2 * (N_ip * get_Qp +dphi*(0.1/(4*sqrt(3)))*randn - y(11)) - 2 * gamma_e * y(15);
dydt(14) = gamma_e^2 * (N_pp * get_Qp - y(10)) - 2 * gamma_e * y(14);  %  dphi multiplies randn because randn gives SD of 1
dydt(15) = gamma_e^2 * (N_ip * get_Qp - y(11)) - 2 * gamma_e * y(15);
dydt(16) = gamma_g^2 * (N_pi * get_Qi - y(12)) - 2 * gamma_g * y(16);
dydt(17) = gamma_g^2 * (N_ii * get_Qi - y(13)) - 2 * gamma_g * y(17);
dydt(18) = (alpha_Na * get_Qp - Na_pump)/tau_Na;
dydt(19) = (g_KNa_0 * (2*y(5))*(1-0.6*y(4))*(1-0.95*y(6)) - y(19))/tau_g;
dydt(20) = (sigma_p_0 - (4*y(4) + 2*y(6)) - y(20))/tau_s;
