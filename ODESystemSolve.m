function [ solution ] = ODESystemSolve( c_sub_i, f_sub_g_of_t, ...
f_sub_m_of_t, I_sub_m_of_t, log_k_sub_12, log_k_sub_a1, log_k_sub_a2, ...
log_k_sub_a3, log_S_sub_t, log_S_sub_d, log_S_sub_e, F_01, ...
EGP_0_minus_F_01, log_k_sub_e, log_k_sub_is1, log_k_sub_is2, ...
log_k_sub_if, log_k_sub_m, log_d, log_Q_sub_b, p_sub_i, p_sub_m, ...
u_sub_i_of_t, G_of_t,V_sub_i,w,V,CHO,t)

% ODESystemSolve is where the system of differential equations is solved 
% for as a system, for now, using ode45 solver

 % c_sub_i, [U/min] background insulin appearance
 % f_sub_g_of_t, [umol/kg/min] additive time-varying piecewise linear flux
 % f_sub_m_of_t, [Unitless], time-varying piecewise linear function
 % I_sub_m_of_t, [Unitless], multiplicative time-varying piecewise linear 
 % function
 % log_k_sub_12, [log(1/min)], natural logarithm of transfer rate constant 
 % from non-accessible to accessible compartment
 % log_k_sub_a1, [log(1/min)], natural logarithm of fractional 
 % deactivation rate constant
 % log_k_sub_a2, [log(1/min)], natural logarithm of fractional 
 % deactivation rate constant
 % log_k_sub_a3, [log(1/min)], natural logarithm of fractional 
 % deactivation rate constant 
 % log_S_sub_t, [log(10^-4/min/mU/L)], natural logarithm of insulin 
 % sensitivity of glucose distribution/transport
 % log_S_sub_d, [log(10^-4/min/mU/L)], natural logarithm of insulin 
 % sensitivity of glucose disposal
 % log_S_sub_e, [log(10^-4/mU/L)], natural logarithm of insulin 
 % sensitivity of endogenous glucose production
 % F_01, [umol/kg/min], noninsulin dependent glucose utilization
 % EGP_0, [mmol/min], endogenous glucose production extrapolated to zero insulin
 % concentration
 % log_k_sub_e, [log(1/min)], natural logarithm of fractional clearance 
 % rate
 % log_k_sub_is1, [log(1/min)], natural logarithm of fractional transfer 
 % rate
 % log_k_sub_is2, [log(1/min)], natural logarithm of fractional transfer 
 % rate
 % log_k_sub_if, [log(1/min)], natural logarithm of shared fractional 
 % transfer rate
 % log_k_sub_m, [log(1/min)], natural logarithm of transfer rate
 % log_d, [log(min)], natural logarithm of delay
 % log_Q_sub_b, [log(U)], natural logarithm of insulin on board due to 
 % preceding insulin delivery
 % p_sub_i, [Unitless], portion of subcutaneous insulin absorbed through the slow
 % channel
 % p_sub_m, [Unitless], portion of meal carbohydrates absorbed through the first
 % channel
 % u_sub_i_of_t, [U/min], insulin infusion rate
 % G_of_t, [mmol/L], plasma blood glucose (presumed to be from regular 
 % blood glucose meter)
 
 % y(1) = Q_sub_is1_of_0, [U], insulin mass in the first compartment in 
 % the slow channel
 y(1) = (((u_sub_i_of_t)*(p_sub_i))/(exp(log_k_sub_is1)));
 
 % y(2) = Q_sub_is2_of_0, [U], insulin mass in the second compartment in
 % the slow channel
 y(2) = (((y(1))*(exp(log_k_sub_is1)))/(exp(log_k_sub_is2)))+ ...
     ((exp(log_Q_sub_b))*p_sub_i);
 
 % y(3) = Q_sub_if1_of_0, [U], insulin mass in the first compartment in
 % the fast channel
 y(3) = (((u_sub_i_of_t)*(1-p_sub_i))/(exp(log_k_sub_if)));
 
 % y(4) = Q_sub_if2_of_0, [U], insulin mass in the second compartment in
 % the fast channel
 y(4) = (y(3)) + ((exp(log_Q_sub_b))*(1-p_sub_i));
 
 % y(5) = Q_sub_i_of_0, [U], insulin mass in the plasma
 y(5) = (((((y(2))*(exp(log_k_sub_is2)))+((y(4))*(exp(log_k_sub_if))))*(I_sub_m_of_t))+c_sub_i)/(exp(log_k_sub_e));
 
 % I_sub_p_of_0, [mU/L], plasma insulin concentration
 I_sub_p_of_0 = ((y(5))/((V_sub_i)*w))*(E+06);
 
 % y(6) = x_sub_1_of_0, [1/min], (remote) effect of insulin on glucose
 % distribution/transport
 y(6)= (exp(log_S_sub_t))*I_sub_p_of_0;
 
 % y(7) = x_sub_2_of_0, [1/min], (remote) effect of insulin on glucose
 % disposal
 y(7) = (exp(log_S_sub_d))*I_sub_p_of_0;
 
 % y(8) = x_sub_3_of_0, [1/min], (remote) effect of insulin on 
 % endogenous glucose production
 y(8) = (exp(log_S_sub_e))*I_sub_p_of_0;
 
 % y(9) = Q_sub_1_of_0, [umol/kg], glucose mass in the accessible
 % (where measurements are made) compartment
 y(9) = (G_of_t)*V;

 % y(10) = Q_sub_2_of_0, [umol/kg], glucose mass in the nonaccessible 
 % compartment
 y(10) = ((y(9))*(y(6)))/((exp(log_k_sub_12))+(y(7)));
 
 % y(11) = p_sub_i, [Unitless], portion of subcutaneous insulin absorbed
 % through the slow channel
 y(11) = p_sub_i;
 
 % y(12) = exp(log_k_sub_is1), [1/min], fractional transfer rate parameter
 y(12) = exp(log_k_sub_is1);
 
 % y(13) = exp(log_k_sub_is2), [1/min], fractional transfer rate parameter
 y(13) = exp(log_k_sub_is2);
 
 % y(14) = exp(log_Q_sub_b), [U], insulin on board due to a preceding
 % insulin delivery
 y(14) = exp(log_Q_sub_b);
 
 % y(15) = exp(log_k_sub_if), [1/min], shared fractional transfer rate
 y(15) = exp(log_k_sub_if);
 
 % y(16) = exp(log_k_sub_e), [1/min], fractional clearance rate
 y(16) = exp(log_k_sub_e);
 
 % y(17) = c_sub_i, [U/min], background insulin appearance
 y(17) = c_sub_i;
 
 % y(18) = I_sub_m_of_t, [Unitless], multiplicative time-varying piecewise
 % linear function
 y(18) = I_sub_m_of_t;
 
 % y(19) = exp(log_k_sub_a1), [1/min], fractional deactivation rate
 % constant
 y(19) = exp(log_k_sub_a1);
 
 % y(20) = exp(log_k_sub_a2), [1/min], fractional deactivation rate
 % constant
 y(20) = exp(log_k_sub_a2);
 
 % y(21) = exp(log_k_sub_a3), [1/min], fractional deactivation rate
 % constant
 y(21) = exp(log_k_sub_a3);
 
 % y(22) = exp(log_S_sub_t), [10^-4*/min/mU/L], insulin sensitivity of 
 % glucose distribution/transport
 y(22) = exp(log_S_sub_t);
 
 % y(23) = exp(log_S_sub_d), [10^-4/min/mU/L], insulin sensitivity of 
 % glucose disposal
 y(23) = exp(log_S_sub_d);
 
 % y(24) = exp(log_S_sub_e), [10^-4/mU/L], insulin sensitivity of 
 % endogenous glucose production
 y(24) = exp(log_S_sub_e);
 
 % y(25) = F_01, [umol/kg/min], noninsulin dependent glucose utilization
 y(25) = F_01;
 
 % y(26) = EGP_0, [mmol/min], endogenous glucose production extrapolated
 % to zero insulin concentration
 y(26) = EGP_0_minus_F_01 + F_01;
 
 % y(27) = f_sub_g_of_t, [umol/kg/min], additive time-varying piecewise
 % linear flux
 y(27) = f_sub_g_of_t;
 
 % y(28) = k_sub_12, [1/min], transfer rate constant from non-accessible
 % to accessible compartment
 y(28) = exp(log_k_sub_12);
 
 % U_sub_m1_of_t, [umol/kg/min], appearance of the first absorption channel
 % for Gut Glucose Absorption Subsystem
 U_sub_m1_of_t = ((CHO*5551)/w)*(p_sub_m)*((exp(log_k_sub_m)).^2)*t*(exp(-exp(log_k_sub_m)*t));
 
 
 % U_sub_m2_of_r, [umol/kg/min], appearance of the second absorption 
 % channel for Gut Glucose Absorption Subsystem
 
 if t>(exp(log_d))
     U_sub_m2_of_t = ((exp(log_k_sub_m)).^2)*(t-d)*((CHO*5551)/w)*(1-p_sub_m)*(exp((-exp(log_k_sub_m))*(t-d)));
     
 else
     U_sub_m2_of_t = 0;
     
 end
 
 % U_sub_m_of_t, [umol/kg/min], appearance of absorption for Gut Glucose 
 % Absorption Subsystem
U_sub_m_of_t = (U_sub_m1_of_t + U_sub_m2_of_t)*(f_sub_m_of_t); 

% Options for solving system of nonlinear differential equations
options = odeset('RelTol', .00001, ...
                 'NonNegative', [1 2 3 4 5 6 7 8 9 10]);
             
% Solution to system of nonlinear differential equations

[T1, Z1] = ode15(@equations, y, u_sub_i_of_t, V_sub_i, w, U_sub_m_of_t, V, options);

function [dydt] = equations(~, y, u_sub_i_of_t, V_sub_i, w, U_sub_m_of_t, V)
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

dydt(1) = ((u_sub_i_of_t)*(y(11)))-((y(1))*(y(12)));

dydt(2) = ((y(1))*(y(12)))-((y(2))*(y(13)));

dydt(3) = ((u_sub_i_of_t)*(1-(y(11))))-((y(3))*(y(15)));

dydt(4) = ((y(3))*(y(15)))-((y(4))*(y(15)));

dydt(5) = ((((y(2))*(y(13)))+((y(3))*(y(15))))*(y(18)))-((y(5))*(y(16)))+(y(17));

dydt(6) = -((y(19))*(y(6)))+((y(19))*(y(22))*((y(5))/(V_sub_i)*(w)));

dydt(7) = -((y(20))*(y(7)))+((y(20))*(y(23))*((y(5))/(V_sub_i)*(w)));

dydt(8) = -((y(21))*(y(8)))+((y(21))*(y(24))*((y(5))/(V_sub_i)*(w)));

if ((y(8))<1)
dydt(9) = -(y(25))*(((y(9))/V)/(1+((y(9))/V)))-((y(6))*(y(9)))+((y(28))*(y(10)))+((y(26))*(1-(y(8))))+(y(27))+U_sub_m_of_t;

else

    dydt(9) = -(y(25))*(((y(9))/V)/(1+((y(9))/V)))-((y(6))*(y(9)))+((y(28))*(y(10)))+(y(27))+U_sub_m_of_t;
    
end

dydt(10) = ((y(6))*(y(9)))-(((y(28))+(y(7)))*(y(10)));

% Parameters are time invariant
dydt(11) = 0;
dydt(12) = 0;
dydt(13) = 0;
dydt(14) = 0;
dydt(15) = 0;
dydt(16) = 0;
dydt(17) = 0;
dydt(18) = 0;
dydt(19) = 0;
dydt(20) = 0;
dydt(21) = 0;
dydt(22) = 0;
dydt(23) = 0;
dydt(24) = 0;
dydt(25) = 0;
dydt(26) = 0;
dydt(27) = 0;
dydt(28) = 0;

 

end
end
