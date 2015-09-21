function [alpha,mse_pred] = get_lamda_tau(sigma, sigma_e, epsilon, lambda)
% this function is to estimate the theoretical mse ,the optimal tuning
% parameter alpha_tru, optimal regularized paramter tau_s

% input : sigma : variance of  measurement noise
%        sigma_e: variance of the error of side information
%        epsilon:  sparsity    
%        lambda 


% output : tau_s; regularization parameter
%          fMSE_tru : theoretical mse
%          alpha_tru: optimal tuning parameter
syms alpha;
syms sigma_stat;
syms z ;
pi = 3.1416;
% lamda = 1e-2;
% tau = 3;
% case 1:  x = 1    0.064
% 1 : z>alpha-1/sigma_stat         % integral region
%  sigma_stat^2*(z-alpha)^2    % the function
% 0 : (-1)*alpha-1/sigma_stat < z < alpha-1/sigma_stat
% 1
% -1 : z < (-1)*alpha-1/sigma_stat
% % sigma_stat^2*(z+alpha)^2
func1_1 =   1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z-alpha)^2 ;
a1 = int(func1_1, z, alpha-1/sigma_stat,inf);
func1_2 =  1/sqrt(2*pi)*exp(-z^2/2);
a2 = int(func1_2, z, (-1)*alpha-1/sigma_stat, alpha-1/sigma_stat);
func1_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
a3 = int(func1_3, z, -inf, (-1)*alpha-1/sigma_stat);
% 
% % case 0:  x = 0    0.872
% % 1 : z>alpha
% % sigma_stat^2*(z-alpha)^2
% % 0: 0
% 
% % -1 : z< (-1)*alpha
% % sigma_stat^2*(z+alpha)^2
func2_1 =   1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z-alpha)^2 ;
b1 = int(func2_1, z, alpha, inf);
func2_2 = 0;
b2 = 0;
func2_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
b3 = int(func2_3, z, -inf, (-1)*alpha);
% 
% % case -1: x = -1    0.064
% % 1: z>alpha+1/sigma_stat
% % sigma_stat^2*(z-alpha)^2
% 
% % 0: (-1)*alpha+1/sigma_stat < z < alpha+1/sigma_stat
% % 1
% % -1: z < (-1)*alpha+1/sigma_stat
% % % sigma_stat^2*(z+alpha)^2
% 
func3_1 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z-alpha)^2 ;
c1 = int(func3_1, z, alpha+1/sigma_stat, inf);
func3_2 =  1/sqrt(2*pi)*exp(-z^2/2);
c2 = int(func3_2, z, (-1)*alpha+1/sigma_stat, alpha+1/sigma_stat);
func3_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
c3 = int(func3_3, z, -inf, (-1)*alpha+1/sigma_stat);

fMSE = (a1+a2+a3)*epsilon/2+(b1+b2+b3)*(1-epsilon)+(c1+c2+c3)*epsilon/2;
% fMSE = epision*((sigma_e^2+sigma^2)*sigma_stat^2-sigma^2*sigma_e^2)/(sigma_e^2-sigma_stat^2);

% detection rate
% 1
func =  1/sqrt(2*pi)*exp(-z^2/2);
a1_dr  = int(func, z, (-1)*alpha-1/sigma_stat, alpha-1/sigma_stat);
a1_conj = 1-vpa(a1_dr);

% 0
a2_dr = int(func, z, (-1)*alpha, alpha);
a2_conj = 1-vpa(a2_dr);

% -1
a3_dr = int(func, z, (-1)*alpha+1/sigma_stat, alpha+1/sigma_stat);
a3_conj = 1-vpa(a3_dr);

% 
equ_dr = 0.064*a1_conj+0.872*a2_conj+0.064*a3_conj;

eq1 = sigma_stat^2-(sigma^2+fMSE/epsilon);
% eq1 = lamda/tau-(sigma_e^2+sigma^2+fMSE/epision)/(sigma^2+fMSE/epision)*alpha*sigma_stat;
eq2 = subs(eq1,alpha,lambda);
S = solve(eq2, sigma_stat);
fmse_num = subs(fMSE,alpha,lambda);
mse_pred = abs(double(subs(fmse_num,abs(double(S)))));
% mse_pred = double(fmse_num_1)*abs(double(S));

 %[mse_empirical]= MSE_EMP(lamda,tau,noise_ratio)