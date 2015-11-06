% This is a TEST of ODESystemSolve.m, with TEST values

% clear, clc, close
clear, clc, close

c_sub_i = 0;

f_sub_g_of_t = 0.8;

f_sub_m_of_t = 0.75;

I_sub_m_of_t = 0.8;

log_k_sub_12 = -2.8;

EGP_0 = 16.1;

log_k_sub_a1 = -5.7;

log_k_sub_a2 = -2.9;

log_k_sub_a3 = -3.7;

log_S_sub_t = 3.7;

log_S_sub_d = 2.3;

log_S_sub_e = 6;

F_01 = 12.1;

log_k_sub_e = -2.63;

log_k_sub_is1 = -3.9120;

log_k_sub_is2 = -3.9120;

log_k_sub_if = -2.708;

log_k_sub_m = -3.9;

log_d = 2.3;

log_Q_sub_b = -0.7;

p_sub_i = 0;

p_sub_m = 0;

u_sub_i_of_t = 1.4/60;

G_of_t = 104/18.1;

V_sub_i = 190;

w = 110;

V = 160;

CHO = 120;

t = 15;

EGP_0_minus_F_01 = EGP_0-F_01;



[T, Z] = ODESystemSolve(c_sub_i, f_sub_g_of_t, f_sub_m_of_t, ...
    I_sub_m_of_t, log_k_sub_12, log_k_sub_a1, log_k_sub_a2, log_k_sub_a3, ...
    log_S_sub_t, log_S_sub_d, log_S_sub_e, F_01, EGP_0_minus_F_01, ...
    log_k_sub_e, log_k_sub_is1, log_k_sub_is2, log_k_sub_if, log_k_sub_m, ...
    log_d, log_Q_sub_b, p_sub_i, p_sub_m, u_sub_i_of_t, G_of_t, V_sub_i, ...
    w, V, CHO, t);
