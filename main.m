% This is the main.m file, where the program is executed from
% Toolboxes required: Statistics Toolbox

% clear, clc, close
clear, clc, close

% format long
format long

% User's weight in kilograms [kg]
% dummy variable for now
w = 110; 
% [kg]

% Infusion rate in units per hour [U/h]
% dummy variable for now
IR_h = 1.4; 
% [U/h]

% Infusion rate in units per min [U/min]
% dummy variable for now
IR_min = 1.4/60; 
% [U/min]

% CHO [g], carbohydrates consumed
% dummy value, for now
CHO = 120; 
% [g]

% c_sub_i, [U/min], background insulin appearance 
% represented as a truncated normal prior distribution for
% Metropolis-Hastings sampling as a Markov chain Monte Carlo method
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers
c_sub_i = mean(truncate(makedist('Normal','mu',0,'sigma',0.04),0,inf)); 
% [U/min]

% log_k_sub_e, [1/min], fractional clearance rate 
% represented as a lognormal prior distribution for Metropolis-Hastings 
% sampling as a Markov chain Monte Carlo method
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers
log_k_sub_e = makedist('Normal','mu',-2.63,'sigma',4.49);
% [1/min]

% log_Q_sub_b, [U], insulin on board due to a preceding insulin delivery
% represented as a lognormal prior distribution for Metropolis-Hastings 
% sampling as a Markov chain Monte Carlo method.
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers
log_Q_sub_b = makedist('Normal','mu',-0.7,'sigma',1.5);
% [U]

% p_sub_i, [unitless], portion of insulin absorbed in the slow channel, in 
% the subcutaneous insulin absorption submodel, represented as a  
% uniform prior distribution for Metropolis-Hastings sampling as a Markov
% chain Monte Carlo method.
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers
p_sub_i = makedist('Uniform',0,1);
% [Unitless]

% p_sub_m, uniform prior distribution
p_sub_m = makedist('Uniform',0,1);

% k_sub_i_mean, contains the mean values of k_sub_is1 (fractional transfer
% rate), k_sub_is2 (fractional transfer rate), k_sub_if (shared fractional
% transfer rate)
k_sub_i_mean = [-3.912,-3.912,-2.708];
% [1/min], [1/min], [1/min]

% omega_sub_i_matrix is the precision matrix for k_sub_is1, k_sub_is2, 
% k_sub_if
omega_sub_i_matrix = [6.25,0,0;0,6.25,0;0,0,6.24];


% Generate random numbers for k_sub_i vector
b = mvnrnd(k_sub_i_mean,inv(omega_sub_i_matrix),100000);

% b_GMM, Gaussian Mixture Model of the multivariate normal randomly
% distributed numbers in the variable b
b_GMM = fitgmdist(b,3,'CovarianceType','diagonal','SharedCov',...
    true);

% Constructing cluster of b_GMM
[idx_b,nlogl_b,P_b,logpdf_b] = cluster(b_GMM,b);

% M_mean, mean for k_sub_m and d
M_mean = [-3.9,2.3];
% [1/min], [min]

% for loop for omega_sub_m_matrix, precision matrix for k_sub_m, and d
if CHO > 0
    % Second value in diagonal for omega_sub_m_matrix
    var = ((100/CHO)^2);
    % [100/CHO] in [100/grams]
    
    % Creating omega_sub_m_matrix
    omega_sub_m_matrix = [12.84,0;0,var];
    
    % Taking inverse of omega_sub_m_matrix
    R1 = inv(omega_sub_m_matrix);
    
    % Generating multivariate normally distributed random numbers
    c = mvnrnd(M_mean,R1,100000);
    
    % c_GMM, Gaussian Mixture Model of the multivariate normal randomly
    % distributed numbers in the variable c
    c_GMM = fitgmdist(c,2,'CovarianceType','diagonal','SharedCov',...
    true);
    
    % Constructing clusters from c_GMM
    [idx_c,nlogl_c,P_c,logpdf_c] = cluster(c_GMM,c);
    
end
% end of for loop

% k_sub_12 data for 6 subject
k_sub_12 = [0.0343,0.0871,0.0863,0.0968,0.0390,0.0458];
% [1/min]

% k_sub_b1 data for 6 subjects, used to find k_sub_a1
k_sub_b1 = [0.0031,0.0157,0.0029,0.0088,0.0007,0.0017];
% [mU/(liter*min^2)]

% S_super_f_sub_IT data for 6 subjects, used to find k_sub_a1
S_super_f_sub_IT = [29.4,18.7,81.2,86.1,72.4,19.1];
% [10^-4/min/mU/L]

% k_sub_a1, deactivation rate constant, data for 6 subjects
k_sub_a1 = (k_sub_b1./S_super_f_sub_IT);
% [1/min]
 

% mean of k_sub_a1 for 6 subjects, deactivation rate constant
mean_k_sub_a1 = mean(k_sub_a1);
% [1/min]

% k_sub_b2, activation rate constant, for 6 subjects, 
% used to find k_sub_a2
k_sub_b2 = [0.0752,0.0231,0.0495,0.0302,0.1631,0.0689];
% [min^-2/mU/L]

% mean of k_sub_b2 for 6 subjects, activation rate constant
mean_k_sub_b2 = mean(k_sub_b2);
% [min^-2/mU/L]

% S_super_f_sub_ID data for 6 subjects, used to find k_sub_a2
S_super_f_sub_ID = [0.9,6.1,20.1,4.7,15.3,2.2];
% [10^-4/min/mU/L]

% k_sub_a2, deactivation rate constant, data for 6 subjects
k_sub_a2 = (k_sub_b2./S_super_f_sub_ID);
% [1/min]

% mean of k_sub_a2, deactivation rate constant, for 6 subjects
mean_k_sub_a2 = mean(k_sub_a2);
% [1/min]

% k_sub_b3, activation rate constant, data for 6 subjects, used to find k_sub_a3
k_sub_b3 = [0.0472,0.0143,0.0691,0.0118,0.0114,0.0285];
% [min^-2/mU/L]

% S_super_f_sub_S_IE data for 6 subjects, used to find k_sub_a3
S_super_f_sub_IE = [401,379,578,720,961,81];
% [10^-4/mU/L]

% k_sub_a3, deactivation rate constant, data for 6 subjects
k_sub_a3 = (k_sub_b3./S_super_f_sub_IE);
% [1/min]

% mean of k_sub_a3, deactivation rate constant, for 6 subjects
mean_k_sub_a3 = mean(k_sub_a3);
% [1/min]

% S_sub_T, insulin sensitivity of insulin distribution/transport,
% data for 6 subjects
S_sub_T = [1.2,10.2,27.1,9.1,16.8,4.7];
% [mL*min^-1*kg^-1/mU/L]

% mean S_sub_T, insulin sensitivity of insulin distribution/transport,
% for 6 subjects
mean_S_sub_T = mean(S_sub_T);
% [mL*min^-1*kg^-1/mU/L]

% S_sub_D, insulin sensitivity of glucose intracellular disposal
% data for 6 subjects
S_sub_D = [4.7,4.7,27.9,15.0,7.3,2.6];
% [mL*min^-1*kg^-1/mU/L]

% mean S_sub_D, insulin sensitivity of glucose intracellular disposal
% for 6 subjects
mean_S_sub_D = mean(S_sub_D);
% [mL*min^-1*kg^-1/mU/L]

% S_sub_E, insulin sensitivity of endogenous glucose production [EGP]
% data for 6 subjects
S_sub_E = [9.0,6.8,14.3,15.2,19.9,1.6];
% [mL*min^-1*kg^-1/mU/L]

% mean S_sub_E, insulin sensitivity of endogenous glucose production [EGP]
mean_S_sub_E = mean(S_sub_E);
% [mL*min^-1*kg^-1/mU/L]

% EGP_0, endogenous glucose production extrapolated to zero insulin
% concentration, data for 6 subjects
EGP_0 = [14.8,14.3,15.6,21.3,20.0,10.5];
% [mmol/min]

% F_01, total non-insulin-dependent glucose flux, data for 6 subjects
F_01 = [12.1,7.5,10.3,11.9,7.1,9.2];
% [mmol/min]

% big_mat, matrix of values for covariance shrinkage
big_mat = [k_sub_12;k_sub_a1;k_sub_a2;k_sub_a3;S_sub_T;S_sub_D;S_sub_E;...
    F_01;EGP_0-F_01];

% P_mean, mean values for big_mat, in log(actual_mean_value) form
% Note: The final value in this vector is incorrect in 
% Stochastic Virtual Population of Subjects With Type 1 Diabetes for the 
% Assessment of Closed-Loop Glucose Controllers and were corrected based
% on data from Partitioning glucose distribution/transport, disposal,
% and endogenous production during IVGTT
% Contacted author, who confirmed error
P_mean = [-2.8,-5.7,-2.9,-3.7,3.7,1.6,6,9.7,6.4];

% cov_shrink, shrinkage of the big_mat data set to a precision matrix
cov_shrink = covshrinkKPM(big_mat',1);

% diag_val, vector of the diagonal values of the covariance matrix
% (inverse of precision matrix)
diag_val = diag((cov_shrink))';

% a, multivariate normal random number generation based on mean P_mean
% and the diagonal of the covariance (inverse precision) matrix
a = mvnrnd(P_mean,diag_val,100000);

% a_GMM, Gaussian Mixture Model of the multivariate normal randomly
% distributed numbers in the variable a
 a_GMM = fitgmdist(a,9,'CovarianceType','diagonal','SharedCov',...
    true);

% Constructing clusters from the a_GMM model
 [idx_a,nlogl_a,P_a,logpdf_a] = cluster(a_GMM,a);

% Concatinating a, b and c matrices
 v = [a,b,c];

% Temporary fit into Gaussian Mixture Model
% THIS WILL NEED TO BE MODIFIED/DELETED
v_GMM = fitgmdist(v,14,'CovarianceType','diagonal','SharedCov',...
  true);

% Creating statistical data of Gaussian Mixture Model of all of the
% parameters
% THIS WILL NEED TO BE MODIFIED/DELETED
 [idx_v,nlogl_v,P_v,logpdf_v] = cluster(v_GMM,v);

% Creating a precision matrix for f_sub_g_of_t that has the dimensions
% equal to the length of random walk for the intrinsic Gaussian Conditional
% AutoRegression, used in Glucose Kinetics Subsystem
f_sub_g_of_t_prec_mat = igmrfprec([15,15],1);

% Creating a precision matrix for f sub m of t that has the dimensions
% equal to the length of random walk for the intrinsic Gaussian Conditional
% AutoRegression 
f_sub_m_of_t_prec_mat = igmrfprec([15,15],1);

% Creating a precision matrix for I sub m of t that has the dimensions
% equal to the length of random walk for the intrinsic Gaussian Conditional
% AutoRegression (intrinsic Gaussian Markov Random Field), used in
% Plasma Insulin Kinetics Subsystem
I_sub_m_of_t_prec_mat = igmrfprec([30,30],1);


% Sigma_sub_g, for f sub g of t, the additive time-varying piecewise flux
% Equivalent to standard deviation (SD) SD = 0.5 [umol/kg/h], used in
% the Glucose Kinetics Subsystem
sigma_sub_g = 16;

% Sigma_sub_i, for f sub m of t, the time-varying piecewise-linear
% function, used in the Gut Glucose Absorption Subsystem
sigma_sub_i_for_f_sub_m_of_t = 100;

% Sigma sub i, for I sub m of t which is the multiplicative time varying 
% piecewise linear function representing diurnal and other time-varying 
% components, used in Plasma Insulin Kinetics Subsystem
sigma_sub_i_for_I_sub_m_of_t = 100;


% V, glucose distribution volume
V = 160;
% [mL/kg]

% V_sub_i, insulin distribution volume
V_sub_i = 190;
% [mL/kg]



