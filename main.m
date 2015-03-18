% This is the main.m file, where the program is executed from
% Toolboxes required: Statistics Toolbox

% clear, clc, close
clear, clc, close

% format long
format long

% CHO [g], carbohydrates consumed
% dummy value, for now
CHO = 120;

% c_sub_i, truncated normal prior distribution
c_sub_i = truncate(makedist('Normal','mu',0,'sigma',(inv(0.04))),0,inf);

% log_k_sub_e, normal prior distribution, with lognormal data as output
log_k_sub_e = makedist('Normal','mu',-2.63,'sigma',(inv(4.49)));

% log_Q_sub_b, normal prior distribution, with lognormal data as output
log_Q_sub_b = makedist('Normal','mu',-0.7,'sigma',(inv(1.5)));

% p_sub_i, uniform prior distribution
p_sub_i = makedist('Uniform',0,1);

% p_sub_m, uniform prior distribution
p_sub_m = makedist('Uniform',0,1);

% k_sub_i_mean, contains the mean values of k_sub_is1, k_sub_is2, k_sub_if
k_sub_i_mean = [-3.912,-3.912,-2.708];

% omega_sub_i_matrix is the precision matrix for k_sub_is1, k_sub_is2, 
% k_sub_if
omega_sub_i_matrix = [6.25,0,0;0,6.25,0;0,0,6.24];

% Generate random numbers for k_sub_i vector
b = mvnrnd(k_sub_i_mean,inv(omega_sub_i_matrix),100000);

% Exponential of random numbers for k_sub_i, for testing correct values
b_exp = exp(b);

% b_GMM, Gaussian Mixture Model of the multivariate normal randomly
% distributed numbers in the variable b
b_GMM = fitgmdist(b,3,'CovarianceType','diagonal','SharedCov',...
    true);

% Constructing cluster of b_GMM
[idx_b,nlogl_b,P_b,logpdf_b] = cluster(b_GMM,b);

% M_mean, mean for k_sub_m and d
M_mean = [-3.9,2.3];

% for loop for omega_sub_m_matrix, precision matrix for k_sub_m, and d
if CHO > 0
    % Second value in diagonal for omega_sub_m_matrix
    var = ((100/CHO)^2);
    
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

% k_sub_b1 data for 6 subjects, used to find k_sub_a1
k_sub_b1 = [0.0031,0.0157,0.0029,0.0088,0.0007,0.0017];

% S_super_f_sub_IT data for 6 subjects, used to find k_sub_a1
S_super_f_sub_IT = [29.4,18.7,81.2,86.1,72.4,19.1];

% k_sub_a1 data for 6 subjects
k_sub_a1 = (k_sub_b1./S_super_f_sub_IT);

% k_sub_b2 data for 6 subjects, used to find k_sub_a2
k_sub_b2 = [0.0752,0.0231,0.0495,0.0302,0.1631,0.0689];

% S_super_f_sub_ID data for 6 subjects, used to find k_sub_a2
S_super_f_sub_ID = [0.9,6.1,20.1,4.7,15.3,2.2];

% k_sub_a2 data for 6 subjects
k_sub_a2 = (k_sub_b2./S_super_f_sub_ID);

% k_sub_b3 data for 6 subjects, used to find k_sub_a3
k_sub_b3 = [0.0472,0.0143,0.0691,0.0118,0.0114,0.0285];

% S_super_f_sub_S_IE data for 6 subjects, used to find k_sub_a3
S_super_f_sub_IE = [401,379,578,720,961,81];

% k_sub_a3 data for 6 subjects
k_sub_a3 = (k_sub_b3./S_super_f_sub_IE);

% S_sub_T data for 6 subjects
S_sub_T = [1.2,10.2,27.1,9.1,16.8,4.7];

% S_sub_D data for 6 subjects
S_sub_D = [4.7,4.7,27.9,15.0,7.3,2.6];

% S_sub_E data for 6 subjects
S_sub_E = [9.0,6.8,14.3,15.2,19.9,1.6];

% EGP_0 data for 6 subjects
EGP_0 = [14.8,14.3,15.6,21.3,20.0,10.5];

% F_01 data for 6 subjects
F_01 = [12.1,7.5,10.3,11.9,7.1,9.2];

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

% prec_shrink, shrinkage of the big_mat data set to a precision matrix
prec_shrink = covshrinkKPM(big_mat',1);

% diag_val, vector of the diagonal values of the covariance matrix
% (inverse of precision matrix)
diag_val = diag(inv(prec_shrink))';

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
v_GMM = fitgmdist(v,14,'CovarianceType','diagonal','SharedCov',...
  true);

% Creating statistical data of Gaussian Mixture Model
 [idx_v,nlogl_v,P_v,logpdf_v] = cluster(v_GMM,v);
  
prec_mat = igmrfprec([10,10],1);

[V,lambda] = eigs(prec_mat);
