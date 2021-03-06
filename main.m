% This is the main.m file, where the program is executed from
% Toolboxes required: Statistics Toolbox, Optimization Toolbox

% clear, clc, close
clear, clc, close

% format long
format long

% Opening data file containing insulin pump, continuous glucose monitor,
% and blood glucose meter readings
fid = fopen('CareLink.csv');

% Separating non-data containing headers and columns from the import
junk = textscan(fid,'%s', 104);

% Getting the data from the headers imported
% M is data matriz
M = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %q %s %s %s %s', 'delimiter', ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;,,,,');

% Indices of the logs cells
indices_cell_string = M{:,1};

% Indices of the logs matrix, in double form
indices_double = str2double(indices_cell_string);

% Date log cells
date_cell_string = M{:,2};

% Time log cells
time_cell_string = M{:,3};

% Timestamp log cells 
timestamp_cell_string = M{:,4};

% New Device Time log cells
new_device_time_cell_string = M{:,5};

% BG reading [mg/dL] cells
BG_reading_cell_string = M{:,6};

% BG reading [mg/dL], in double form
BG_reading_double = str2double(BG_reading_cell_string);

% Linked BG Meter ID cells
linked_BG_meter_ID_cell_string = M{:,7};

% Temp basal amount [U/h] cells
temp_basal_amount_cell_string = M{:,8};

% Temp basal amount [U/h], in double form
temp_basal_amount_double = str2double(temp_basal_amount_cell_string);

% Temp basal type cells
temp_basal_type_cell_string = M{:,9};

% Temp basal duration (HH:MM:SS) cells
temp_basal_duration_cell_string = M{:,10};

% Bolus type cells
bolus_type_cell_string = M{:,11};

% Bolus volume selected [U] cells
bolus_volume_selected_cell_string = M{:,12};

% Bolus volume selected [U], in double form
bolus_volume_selected_double = str2double(bolus_volume_selected_cell_string);

% Bolus volume delivered [U] cells
bolus_volume_delivered_cell_string = M{:,13};

% Bolus volume delivered [U], in double form
bolus_volume_delivered_double = str2double(bolus_volume_delivered_cell_string);

% Programmed bolus duration [HH:MM:SS] cells
programmed_bolus_duration_cells = M{:,14};

% Prime type cells
prime_type_cell_string = M{:,15};

% Prime volume delivered [U] cells
prime_volume_delivered_cell_string = M{:,16};

% Prime volume delivered [U], double form
prime_volume_delivered_double = str2double(prime_volume_delivered_cell_string);

% Suspend cells
suspend_cell_string = M{:,17};

% Rewind cells
rewind_cell_string = M{:,18};

% Bolus wizard estimate [U] cells
bolus_wizard_estimate_cell_string = M{:,19};

% Bolus wizard estimate [U], double form
bolus_wizard_estimate_double = str2double(bolus_wizard_estimate_cell_string);

% Bolus wizard target high BG [mg/dL] cells
bolus_wizard_target_high_BG_cell_string = M{:,20};

% Bolus wizard target high BG [mg/dL], double form
bolus_wizard_target_high_BG_double = str2double(bolus_wizard_target_high_BG_cell_string);

% Bolus wizard target low BG [mg/dL] cells
bolus_wizard_target_low_BG_cell_string = M{:,21};

% Bolus wizard target low BG [mg/dL], double form
bolus_wizard_target_low_BG_double = str2double(bolus_wizard_target_low_BG_cell_string);

% Bolus wizard carb ratio [grams] cells
bolus_wizard_carb_ratio_cell_string = M{:,22};

% Bolus wizard carb ratio [grams], double form
bolus_wizard_carb_ratio_double = str2double(bolus_wizard_carb_ratio_cell_string);

% Bolus wizard insulin sensitivity [mg/dL] cells
bolus_wizard_insulin_sensitivity_cell_string = M{:,23};

% Bolus wizard insulin sensitivity [mg/dL], double form
bolus_wizard_insulin_sensitivity_double = str2double(bolus_wizard_insulin_sensitivity_cell_string);

% Bolus wizard carb input [grams] cells
bolus_wizard_carb_input_cell_string = M{:,24};

% Bolus wizard carb input [grams], double form
bolus_wizard_carb_input_double = str2double(bolus_wizard_carb_input_cell_string);

% Bolus wizard BG input [mg/dL] cells
bolus_wizard_BG_input_cell_string = M{:,25};

% Bolus wizard BG input [mg/dL], double form
bolus_wizard_BG_input_double = str2double(bolus_wizard_BG_input_cell_string);

% Bolus wizard correction estimate [U] cells
bolus_wizard_correction_estimate_cell_string = M{:,26};

% Bolus wizard correction estimate [U], double form
bolus_wizard_correction_estimate_double = str2double(bolus_wizard_correction_estimate_cell_string);

% Bolus wizard food estimate [U] cells
bolus_wizard_food_estimate_cell_string = M{:,27};

% Bolus wizard food estimate [U], double form
bolus_wizard_food_estimate_double = str2double(bolus_wizard_food_estimate_cell_string);

% Bolus wizard active insulin [U] cells
bolus_wizard_active_insulin_cell_string = M{:,28};

% Bolus wizard active insulin [U], double form
bolus_wizard_active_insulin_double = str2double(bolus_wizard_active_insulin_cell_string);

% Alarm cells
Alarm = M{:,29};

% Sensor calibration BG [mg/dL] cells
sensor_calibration_BG_cell_string = M{:,30};

% Sensor calibration BG [mg/dL], double form
sensor_calibration_BG_double = str2double(sensor_calibration_BG_cell_string);

% Sensor glucose [mg/dL] cells
sensor_glucose_cell_string = M{:,31};

% Sensor glucose [mg/dL], double form
sensor_glucose_double = str2double(sensor_glucose_cell_string);

% ISIG value cells
ISIG_value_cell_string = M{:,32};

% ISIG value, double form
ISIG_value_double = str2double(ISIG_value_cell_string);

% Daily insulin total [U] cells
daily_insulin_total_cell_string = M{:,33};

% Daily insulin total [U], double form
daily_insulin_total_double = str2double(daily_insulin_total_cell_string);

% Raw type cells
raw_type_cell_string = M{:,34};

% Raw value cells
raw_value_cell_string = M{:,35};

% Raw ID cells
raw_ID_cell_string = M{:,36};

% Raw upload ID cells
raw_upload_ID_cell_string = M{:,37};

% Raw seq number cells
raw_seq_number_cell_string = M{:,38};

% Raw device type cells
raw_device_type_cell_string = M{:,39};

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
c_sub_i = truncate(makedist('Normal','mu',0,'sigma',0.04),0,inf); 
% [U/min]

% log_k_sub_e, [1/min], fractional clearance rate 
% represented as a lognormal prior distribution for Metropolis-Hastings 
% sampling as a Markov chain Monte Carlo method
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers
log_k_sub_e = makedist('Normal','mu',-2.63,'sigma',4.49); % [log of 1/min]

% log_Q_sub_b, [U], insulin on board due to a preceding insulin delivery
% represented as a lognormal prior distribution for Metropolis-Hastings 
% sampling as a Markov chain Monte Carlo method.
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527,
% right column
log_Q_sub_b = makedist('Normal','mu',-0.7,'sigma',1.5); % [log of U]

% p_sub_i, [unitless], portion of insulin absorbed in the slow channel, in 
% the subcutaneous insulin absorption submodel, represented as a  
% uniform prior distribution for Metropolis-Hastings sampling as a Markov
% chain Monte Carlo method.
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527
% right column
p_sub_i = makedist('Uniform',0,1); % [Unitless]

% p_sub_m, uniform prior distribution
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527,
% right column
p_sub_m = makedist('Uniform',0,1); % [Unitless]

% k_sub_i_mean, contains the population mean values of k_sub_is1 
% (fractional transfer rate), k_sub_is2 (fractional transfer rate), 
% k_sub_if (shared fractional transfer rate)
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527
% middle column
k_sub_i_mean = [-3.912,-3.912,-2.708]; % [1/min], [1/min], [1/min]

% omega_sub_i_matrix is the precision matrix, which is inverse of the 
% covariance matrix, for k_sub_is1, k_sub_is2, k_sub_if
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527
% middle column
omega_sub_i_matrix = [6.25,0,0;0,6.25,0;0,0,6.24];

% Generate random numbers for k_sub_i vector
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers
b = mvnrnd(k_sub_i_mean,inv(omega_sub_i_matrix),100000);

% b_GMM, Gaussian Mixture Model of the multivariate normal randomly
% distributed numbers in the variable b
b_GMM = fitgmdist(b,3,'CovarianceType','diagonal','SharedCov',...
    true);

% Constructing cluster of b_GMM
[idx_b,nlogl_b,P_b,logpdf_b] = cluster(b_GMM,b);

% M_mean, population mean values for k_sub_m and d
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527
% middle column
M_mean = [-3.9,2.3]; % [1/min], [min]

% for loop for omega_sub_m_matrix, precision matrix for k_sub_m, and d
if CHO > 0
    % Second value in diagonal for omega_sub_m_matrix
    % From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
    % for the Assessment of Closed Loop Glucose Controllers, page 3527,
    % middle column
    var = ((100/CHO)^2); % [100/CHO] in [100*grams^-1]
    
    % Creating omega_sub_m_matrix
    % From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
    % for the Assessment of Closed Loop Glucose Controllers, page 3527,
    % middle column
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

% k_sub_12, transfer rate constant from non-accessible to accessible 
% compartment, data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 3
k_sub_12 = [0.0343,0.0871,0.0863,0.0968,0.0390,0.0458]; % [1/min]

% k_sub_b1, activation rate constant, data for 6 subjects, used to find 
% k_sub_a1
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 4
k_sub_b1 = [0.0031,0.0157,0.0029,0.0088,0.0007,0.0017]; % [(min^-2)/mU/L]

% S_super_f_sub_IT data for 6 subjects, used to find k_sub_a1
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 7
S_super_f_sub_IT = [29.4,18.7,81.2,86.1,72.4,19.1]; % [10^-4/min/mU/L]

% k_sub_a1, deactivation rate constant, data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, equation used is from description 
% underneath Table 1 on page E998
k_sub_a1 = (k_sub_b1./S_super_f_sub_IT); % [1/min]

% k_sub_b2, activation rate constant, for 6 subjects, 
% used to find k_sub_a2
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 5
k_sub_b2 = [0.0752,0.0231,0.0495,0.0302,0.1631,0.0689]; % [min^-2/mU/L]

% S_super_f_sub_ID data for 6 subjects, used to find k_sub_a2
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 8
S_super_f_sub_ID = [0.9,6.1,20.1,4.7,15.3,2.2]; % [10^-4/min/mU/L]

% k_sub_a2, deactivation rate constant, data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, equation used is from description 
% underneath Table 1 on page E998
k_sub_a2 = (k_sub_b2./S_super_f_sub_ID); % [1/min]

% k_sub_b3, activation rate constant, data for 6 subjects, used to find 
% k_sub_a3
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 6
k_sub_b3 = [0.0472,0.0143,0.0691,0.0118,0.0114,0.0285]; % [min^-2/mU/L]

% S_super_f_sub_S_IE data for 6 subjects, used to find k_sub_a3
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 9
S_super_f_sub_IE = [401,379,578,720,961,81]; % [10^-4/mU/L]

% k_sub_a3, deactivation rate constant, data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, equation used is from description
% underneath Table 1 on page E998
k_sub_a3 = (k_sub_b3./S_super_f_sub_IE); % [1/min]

% S_sub_T, insulin sensitivity of insulin distribution/transport,
% data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E999, Table 2, column 2
S_sub_T = [1.2,10.2,27.1,9.1,16.8,4.7]; % [mL*min^-1*kg^-1/mU/L]

% S_sub_D, insulin sensitivity of glucose intracellular disposal
% data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E999, Table 2, column 3
S_sub_D = [4.7,4.7,27.9,15.0,7.3,2.6]; % [mL*min^-1*kg^-1/mU/L]

% S_sub_E, insulin sensitivity of endogenous glucose production [EGP]
% data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E999, Table 2, column 4
S_sub_E = [9.0,6.8,14.3,15.2,19.9,1.6]; % [mL*min^-1*kg^-1/mU/L]

% EGP_0, endogenous glucose production extrapolated to zero insulin
% concentration, data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 12
EGP_0 = [14.8,14.3,15.6,21.3,20.0,10.5]; % [mmol/min]

% F_01, total non-insulin-dependent glucose flux, data for 6 subjects
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, page E998, Table 1, column 10
F_01 = [12.1,7.5,10.3,11.9,7.1,9.2]; % [mmol/min]

% big_mat, matrix of values for covariance shrinkage
% From: Partitioning glucose distribution/transport, disposal, and 
% endogenous production during IVGTT, pages E998 (Table 1) 
% and E999 (Table 2)
big_mat = [k_sub_12;k_sub_a1;k_sub_a2;k_sub_a3;S_sub_T;S_sub_D;S_sub_E;...
    F_01;EGP_0-F_01];

% P_mean, mean values for big_mat, in log(actual_mean_value) form
% From: Stochastic Virtual Population of Subjects With
% Type 1 Diabetes for the Assessment of Closed-Loop
% Glucose Controllers, page 3527
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

% Sigma_sub_g, for f sub g of t, the additive time-varying piecewise flux
% Equivalent to standard deviation (SD) SD = 0.5 [umol/kg/h], used in
% the Glucose Kinetics Subsystem
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527
sigma_sub_g = 16;

% Sigma_sub_i, for f sub m of t, the time-varying piecewise-linear
% function, used in the Gut Glucose Absorption Subsystem
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527
sigma_sub_i_for_f_sub_m_of_t = 100;

% Sigma sub i, for I sub m of t which is the multiplicative time varying 
% piecewise linear function representing diurnal and other time-varying 
% components, used in Plasma Insulin Kinetics Subsystem
% From: Stochastic Virtual Population of Subjects with Type 1 Diabetes
% for the Assessment of Closed Loop Glucose Controllers, page 3527
sigma_sub_i_for_I_sub_m_of_t = 100;

% V, glucose distribution volume
% From: Stochastic Virtual Population of Subjects With
% Type 1 Diabetes for the Assessment of Closed-Loop
% Glucose Controllers, page 3526
V = 160; % [mL/kg]

% V_sub_i, insulin distribution volume
% From: Stochastic Virtual Population of Subjects With
% Type 1 Diabetes for the Assessment of Closed-Loop
% Glucose Controllers, page 3525
V_sub_i = 190; % [mL/kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Initializing parameter guesses                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial guesses for c_sub_i, insulin appearance
c_sub_i_initial_guess_1 = 0.0; % [U/min]
c_sub_i_initial_guess_2 = 0.0; % [U/min]
c_sub_i_initial_guess_3 = 0.0; % [U/min]

% Initial guesses for f_sub_g_of_t, additive time-varying piecewise linear 
% flux AKA "glucose flux"
% Consider having an initial guess for each random walk
f_sub_g_of_t_initial_guess_1 = 0; % [umol/kg/min]
f_sub_g_of_t_initial_guess_2 = 0; % [umol/kg/min]
f_sub_g_of_t_initial_guess_3 = 0; % [umol/kg/min]

% Initial guesses of dinnertime f_sub_m_of_t, time-varying piecewise linear 
% function
% Consider having an initial guess for each random walk
% Link this with the inital dinnertime meal bolus
f_sub_m_of_t_dinnertime_initial_guess_1 = 0; % [Unitless]
f_sub_m_of_t_dinnertime_initial_guess_2 = 0; % [Unitless]
f_sub_m_of_t_dinnertime_initial_guess_3 = 0; % [Unitless]

% Initial guesses of breakfast f_sub_m_of_t, time-varying piecewise linear
% function
% Consider having an initial guess for each random walk
% Link this with the initial breakfast meal bolus
f_sub_m_of_t_breakfast_initial_guess_1 = 0; % [Unitless]
f_sub_m_of_t_breakfast_initial_guess_2 = 0; % [Unitless]
f_sub_m_of_t_breakfast_initial_guess_3 = 0; % [Unitless]


% Initial guesses of lunchtime f_sub_m_of_t, time-varying piecewise linear
% function
% Consider having an initial guess for each random walk
% Link this with the initial lunchtime meal bolus
f_sub_m_of_t_lunchtime_initial_guess_1 = 0; % [Unitless]
f_sub_m_of_t_lunchtime_initial_guess_2 = 0; % [Unitless]
f_sub_m_of_t_lunchtime_initial_guess_3 = 0; % [Unitless]

% Initial guesses for I_sub_m_of_t, multiplicative
% time-varying piecewise-linear function describing diurnal and other
% time-varying components of insulin-kinetics
% Consider having an initial guess for each random walk
I_sub_m_of_t_initial_guess_1 = 0; % [Unitless]
I_sub_m_of_t_initial_guess_2 = 0; % [Unitless]
I_sub_m_of_t_initial_guess_3 = 0; % [Unitless]

% Initial guesses for log_k_sub_12, the natural logarithm of the transfer
% rate constant from the non-accessible to accessible compartment
% Initial guess equals the population mean value of the population 
% distribution for the parameter
log_k_sub_12_initial_guess_1 = -2.8; % [log(1/min)]
log_k_sub_12_initial_guess_2 = -2.8; % [log(1/min)]
log_k_sub_12_initial_guess_3 = -2.8; % [log(1/min)]

% Initial guesses for log_k_sub_a1, the natural logarithm of a fractional 
% deactivation rate constant
% Initial guess equals the population mean value of the population 
% distribution for the parameter
log_k_sub_a1_initial_guess_1 = -5.7; % [log(1/min)]
log_k_sub_a1_initial_guess_2 = -5.7; % [log(1/min)]
log_k_sub_a1_initial_guess_3 = -5.7; % [log(1/min)]

% Initial guesses for log_k_sub_a2, the natural logarithm of a fractional 
% deactivation rate constant
% Initial guess equals the population mean value of the population 
% distribution for the parameter
log_k_sub_a2_initial_guess_1 = -2.9; % [log(1/min)]
log_k_sub_a2_initial_guess_2 = -2.9; % [log(1/min)]
log_k_sub_a2_initial_guess_3 = -2.9; % [log(1/min)]

% Initial guesses for log_k_sub_a3, the natural logarithm of a fractional 
% deactivation rate constant
% Initial guess equals the population mean value of the population 
% distribution for the parameter
log_k_sub_a3_initial_guess_1 = -3.7; % [log(1/min)]
log_k_sub_a3_initial_guess_2 = -3.7; % [log(1/min)]
log_k_sub_a3_initial_guess_3 = -3.7; % [log(1/min)]

% Initial guesses for log_S_sub_t, the natural logarithm of the insulin 
% sensitivity of insulin distribution/transport
% Initial guess equals the population mean value of the population 
% distribution for the parameter
log_S_sub_t_initial_guess_1 = 3.7; % [log((10^-4)*/min/mU/L)]
log_S_sub_t_initial_guess_2 = 3.7; % [log((10^-4)*/min/mU/L)]
log_S_sub_t_initial_guess_3 = 3.7; % [log((10^-4)*/min/mU/L)]

% Initial guesses for log_S_sub_d, the natural logarithm of the insulin 
% sensitivity of glucose disposal
% Initial guess equals the population mean value of the population 
% distribution for the parameter
log_S_sub_d_initial_guess_1 = 1.6; % [log((10^-4)*/min/mU/L)]
log_S_sub_d_initial_guess_2 = 1.6; % [log((10^-4)*/min/mU/L)]
log_S_sub_d_initial_guess_3 = 1.6; % [log((10^-4)*/min/mU/L)]

% Initial guesses for log_S_sub_e, the natural logarithm of the insulin 
% sensitivity of endogenous glucose production (EGP)
log_S_sub_e_initial_guess_1 = 6; % [log((10^-4)*/mU/L)]
log_S_sub_e_initial_guess_2 = 6; % [log((10^-4)*/mU/L)]
log_S_sub_e_initial_guess_3 = 6; % [log((10^-4)*/mU/L)]

% Initial guesses for F_sub_01, the noninsulin dependent glucose utilization
F_sub_01_initial_guess_1 = 7.3; % [umol*(kg^-1)*(min^-1)]
F_sub_01_initial_guess_2 = 7.3; % [umol*(kg^-1)*(min^-1)]
F_sub_01_initial_guess_3 = 7.3; % [umol*(kg^-1)*(min^-1)]

% Initial guesses for EGP_sub_0, endogenous insulin production extrapolated
% to zero insulin concentration
EGP_sub_0_initial_guess_1 = 26.3; % [umol*(kg^-1)*(min^-1)]
EGP_sub_0_initial_guess_2 = 26.3; % [umol*(kg^-1)*(min^-1)]
EGP_sub_0_initial_guess_3 = 26.3; % [umol*(kg^-1)*(min^-1)]

% Initial guesses for log_k_sub_e, the natural logarithm of the fractional 
% clearance rate of plasma insulin
log_k_sub_e_initial_guess_1 = -2.63; % [log(1/min)]
log_k_sub_e_initial_guess_2 = -2.63; % [log(1/min)]
log_k_sub_e_initial_guess_3 = -2.63; % [log(1/min)]

% Initial guesses for log_k_sub_is1, the natural logarithm of the fractional 
% transfer rate parameter
log_k_sub_is1_initial_guess_1 = -3.912; % [log(1/min)]
log_k_sub_is1_initial_guess_2 = -3.912; % [log(1/min)]
log_k_sub_is1_initial_guess_3 = -3.912; % [log(1/min)]

% Initial guesses for log_k_sub_is2, the natural logarithm of the fractional 
% transfer rate parameter
log_k_sub_is2_initial_guess_1 = -3.912; % [log(1/min)]
log_k_sub_is2_initial_guess_2 = -3.912; % [log(1/min)]
log_k_sub_is2_initial_guess_3 = -3.912; % [log(1/min)]

% Initial guesses for log_k_sub_if, the natural logarithm of the shared 
% fractional transfer rate parameter
log_k_sub_if_initial_guess_1 = -2.708; % [log(1/min)]
log_k_sub_if_initial_guess_2 = -2.708; % [log(1/min)]
log_k_sub_if_initial_guess_3 = -2.708; % [log(1/min)]

% Initial guesses for log_k_sub_m, the natural logarithm of the transfer rate 
% parameter
log_k_sub_m_initial_guess_1 = -3.9; % [log(1/min)]
log_k_sub_m_initial_guess_2 = -3.9; % [log(1/min)]
log_k_sub_m_initial_guess_3 = -3.9; % [log(1/min)]

% Initial guesses for log_of_d, the natural logarithm of the delay associated 
% with the second absorption channel
log_sub_of_d_initial_guess_1 = 2.3; % [log(min)]
log_sub_of_d_initial_guess_2 = 2.3; % [log(min)]
log_sub_of_d_initial_guess_3 = 2.3; % [log(min)]

% Initial guesses for log_Q_sub_b, the natural logarithm of the insulin on 
% board due to a preceding insulin delivery
log_Q_sub_b_initial_guess_1 = -0.7; % [log(U)]
log_Q_sub_b_initial_guess_2 = -0.7; % [log(U)]
log_Q_sub_b_initial_guess_3 = -0.7; % [log(U)]

% Initial guesses for log_p_sub_i, the natural logarithm for the portion of 
% subcutaneous insulin absorbed through the slow channel
log_p_sub_i_initial_guess_1 = 0; % [log(Unitless)] 
log_p_sub_i_initial_guess_2 = 0; % [log(Unitless)]
log_p_sub_i_initial_guess_3 = 0; % [log(Unitless)]

% Initial guesses for log_p_sub_m, the natural logarithm for the portion of 
% meal carbohydrates absorbed in the first channel
 log_p_sub_m_initial_guess_1 = 0; % [log(Unitless)]
 log_p_sub_m_initial_guess_2 = 0; % [log(Unitless)]
 log_p_sub_m_initial_guess_3 = 0; % [log(Unitless)]
 




