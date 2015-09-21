function [y] = generate_A(z,chosen_ind,m,n)

N = m*n;

z = reshape(z,m,n);

alpha = idct2(z);

y_vec = alpha(:);

% indSample = randperm(N);

% chosen_ind = indSample(1:M);

y = y_vec(chosen_ind);