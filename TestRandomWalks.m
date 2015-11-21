% Random walk test

% clear, clc, close
clear, clc, close

f_sub_g_of_t_precision_matrix = igmrfprec([4,4],1);

f_sub_g_of_t_dist = makedist('Normal','mu',0,'sigma',16);

x_sub_i_for_f_sub_g_of_t = random(f_sub_g_of_t_dist,[16,1]);

[U,S,V] = svds(f_sub_g_of_t_precision_matrix);

rank_f_sub_g_of_t_precision_matrix = rank(U*S*V');



