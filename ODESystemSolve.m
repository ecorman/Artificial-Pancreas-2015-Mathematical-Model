function [ solution ] = ODESystemSolve( c_sub_i, f_sub_g_of_t, ...
f_sub_m_of_t, I_sub_m_of_t, log_k_sub_12, log_k_sub_a1, log_k_sub_a2, ...
log_k_sub_a3, log_S_sub_t, log_S_sub_d, log_S_sub_e, F_01, ...
EGP_0_minus_F_01, log_k_sub_e, log_k_sub_is1, log_k_sub_is2, ...
log_k_sub_if, log_k_sub_m, log_d, log_Q_sub_b, p_sub_i, p_sub_m, ...
u_sub_i_of_t, G_of_t)

% ODESystemSolve is where the system of differential equations is solved 
% for as a system, for now, using ode45 solver

 % y(1) = c_sub_i, [U/min] background insulin appearance
 % y(2) = f_sub_g_of_t, [umol/kg/min] additive time-varying piecewise linear flux
 % y(3) = f_sub_m_of_t, [Unitless], time-varying piecewise linear function
 % y(4) = I_sub_m_of_t, [Unitless], multiplicative time-varying piecewise linear 
 % function
 % y(5) = log_k_sub_12, [log(1/min)], natural logarithm of transfer rate constant 
 % from non-accessible to accessible compartment
 % y(6) = log_k_sub_a1, [log(1/min)], natural logarithm of fractional 
 % deactivation rate constant
 % y(7) = log_k_sub_a2, [log(1/min)], natural logarithm of fractional 
 % deactivation rate constant
 % y(8) = log_k_sub_a3, [log(1/min)], natural logarithm of fractional 
 % deactivation rate constant 
 % y(9) = log_S_sub_t, [log(10^-4/min/mU/L)], natural logarithm of insulin 
 % sensitivity of glucose distribution/transport
 % y(10) = log_S_sub_d, [log(10^-4/min/mU/L)], natural logarithm of insulin 
 % sensitivity of glucose disposal
 % y(11) = log_S_sub_e, [log(10^-4/mU/L)], natural logarithm of insulin 
 % sensitivity of endogenous glucose production
 % y(12) = F_01, [umol/kg/min], noninsulin dependent glucose utilization
 % y(13) = EGP_0, [mmol/min], endogenous glucose production extrapolated to zero insulin
 % concentration
 % y(14) = log_k_sub_e, [log(1/min)], natural logarithm of fractional clearance 
 % rate
 % y(15) = log_k_sub_is1, [log(1/min)], natural logarithm of fractional transfer 
 % rate
 % y(16) = log_k_sub_is2, [log(1/min)], natural logarithm of fractional transfer 
 % rate
 % y(17) = log_k_sub_if, [log(1/min)], natural logarithm of shared fractional 
 % transfer rate
 % y(18) = log_k_sub_m, [log(1/min)], natural logarithm of transfer rate
 % y(19) = log_d, [log(min)], natural logarithm of delay
 % y(20) = log_Q_sub_b, [log(U)], natural logarithm of insulin on board due to 
 % preceding insulin delivery
 % y(21) = p_sub_i, [Unitless], portion of subcutaneous insulin absorbed through the slow
 % channel
 % y(22) = p_sub_m, [Unitless], portion of meal carbohydrates absorbed through the first
 % channel
 % y(23) = u_sub_i_of_t, [U/min], insulin infusion rate
 % y(24) = G_of_t, [mmol/L], plasma blood glucose (presumed to be from regular 
 % blood glucose meter)

end

function [dydt] = equations(~,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subcutaneous Insulin Absorption Subsystem
% 
% Equations to be solved for:
%
% dydt(1) = Q_dot_sub_is1_of_t
% dydt(2) = Q_dot_sub_is2_of_t
%
% Q_sub_is1_of_t, [U], is the insulin mass in the first compartment of
% the slow channel
% Q_sub_is2_of_t, [U], is the insulin mass in the second compartment of
% the slow channel
%
% dydt(3) = Q_dot_sub_if1_of_t
% dydt(4) = Q_dot_sub_if2_of_t
%
% Q_sub_if1_of_t, [U], is the insulin mass in the first compartment of
% the fast channel
% Q_sub_if2_of_t, [U], is the insulin mass in the second compartment of
% the fast channel
%
% Initial Conditions
%
% Q_sub_is1_of_0 = u_sub_i_of_0*p_sub_i/k_sub_is1
% Q_sub_is2_of_0 = Q_is1_of_0*k_sub_is1
%
% Input values
% 
% u_sub_i_of_t, [U/min], insulin infusion rate
%
% Time-Invariant Parameters
%
% k_sub_is1, [1/min], fractional deactivation rate
% k_sub_is2, [1/min], fractional deactivation rate
% p_sub_i, [Unitless], portion of subcutaneous insulin absorbed through
% the slow channel
% Q_sub_b, [U], insulin on board due to preceding insulin delivery
% k_sub_if, [1/min], shared fractional transfer rate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plasma (Blood) Insulin Kinetics Subsystem
%
% Equation to be solved for:
%
% dydt(5) = Q_dot_sub_i_of_t
%
% Q_sub_i_of_t, [U], is the insulin mass in the plasma (blood)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Insulin Action Subsystem
%
% Equations to be solved for:
%
% dydt(6) = x_dot_sub_1_of_t
% dydt(7) = x_dot_sub_2_of_t
% dydt(8) = x_dot_sub_3_of_t
%
% x_sub_1_of_t, [1/min], represents the (remote) effect of insulin on  
% glucose distribution/transport
% x_sub_2_of_t, [1/min], represents the (remote) effect of insulin on
% glucose disposal
% x_sub_3_of_t, [1/min], represents the (remote) effect of insulin on
% the endogenous insulin production
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Glucose Kinetics Subsystem
%
% Equations to be solved for:
%
% dydt(9) = Q_dot_sub_1_of_t
% dydt(10) = Q_dot_sub_2_of_t
%
% Q_sub_1_of_t, [umol/kg/min], glucose mass in the accessible (where 
% measurements are made) compartment
% Q_sub_2_of_t, [umol/kg/min], glucose mass in the nonaccessible compartment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dydt = y;





end
