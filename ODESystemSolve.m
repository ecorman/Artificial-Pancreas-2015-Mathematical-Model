function [ solution ] = ODESystemSolve( c_sub_i, f_sub_g_of_t, ...
f_sub_m_of_t, I_sub_m_of_t, log_k_sub_12, log_k_sub_a1, log_k_sub_a2, ...
log_k_sub_a3, log_S_sub_t, log_S_sub_d, log_S_sub_e, F_01, ...
EGP_0_minus_F_01, log_k_sub_e, log_k_sub_is1, log_k_sub_is2, ...
log_k_sub_if, log_k_sub_m, log_d, log_Q_sub_b, p_sub_i, p_sub_m, ...
u_of_t, G_of_t)

% ODESystemSolve is where the system of differential equations is solved 
% for as a system, for now, using ode45 solver

 % c_sub_i, background insulin appearance, [U/min]
 
 % f_sub_g_of_t, additive time-varying piecewise linear flux, [umol/kg/min]
 
 % f_sub_m_of_t, time-varying piecewise linear function, [Unitless]
 
 % I_sub_m_of_t, multiplicative time-varying piecewise linear function,
 % [Unitless]
 
 % log_k_sub_12, log of transfer rate constant from non-accessible to 
 % accessible compartment, [log of 1/min]
 
 % log_k_sub_a1, log of fractional deactivation rate constant, 
 % [log of 1/min]
 
 % log_k_sub_a2, log of fractional deactivation rate constant,
 % [log of 1/min]
 
 % log_k_sub_a3, log of fractional deactivation rate constant, 
 % [log of 1/min]
 
 % log_S_sub_t, log of insulin sensitivity of glucose 
 % distribution/transport, [log of 10^-4/min/mU/L]
 
 % log_S_sub_d, log of insulin sensitivity of glucose disposal, 
 % [log of 10^-4/min/mU/L]
 
 % log_S_sub_e, log of insulin sensitivity of endogenous glucose 
 % production, [log of 10^-4/mU/L]
 
 % F_01, noninsulin dependent glucose utilization, [umol/kg/min]
 
 % EGP_0, endogenous glucose production extrapolated to zero insulin
 % concentration, [mmol/min]
 
 % log_k_sub_e, log of fractional clearance rate, [log of 1/min]
 
 % log_k_sub_is1, log of fractional transfer rate, [log of 1/min]
 
 % log_k_sub_is2, log of fractional transfer rate, [log of 1/min]
 
 % log_k_sub_if, log of shared fractional transfer rate, [log of 1/min]
 
 % log_k_sub_m, log of transfer rate, [log of 1/min]
 
 % log_d, log of delay, [log of min]
 
 % log_Q_sub_b, log of insulin on board due to preceding insulin delivery,
 % [log of Units]
 
 % p_sub_i, portion of subcutaneous insulin absorbed through the slow
 % channel, [Unitless]
 
 % p_sub_m, portion of meal carbohydrates absorbed through the first
 % channel, [Unitless]




end

