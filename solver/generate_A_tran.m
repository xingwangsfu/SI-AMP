function [x] = generate_A_tran(z, chosen_ind, m,n)

N = m*n;

y_vec = zeros(N,1);

y_vec(chosen_ind) = z;

y_matrix = reshape(y_vec,m,n);

x_matrix = dct2(y_matrix);

x = x_matrix(:);