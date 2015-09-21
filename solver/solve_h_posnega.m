function [h_scalar,fMSE_num] = solve_h_posnega(epsilon, mse_object, alpha)

% this function is to find the c-nearly-least-favorable distribution.
% input: epsilon: sparsity
%         mse_object: objective mse
%         alpha: threshold, tuning parameter
%        sigma_stat: variance of measurement noise

% output: h_scalar: the position where the mass places
% syms z ;
% syms h_posnega;
% % pi = 3.1416;
% % lamda = 1e-2;
% % tau = 3;
% % case 1:  x = h_posnega    epision/2
% % 1 : z>alpha-h_posnega/sigma_stat         % integral region
% %  sigma_stat^2*(z-alpha)^2    % the function
% % 0 : (-1)*alpha-h_posnega/sigma_stat < z < alpha-h_posnega/sigma_stat
% % 1
% % -1 : z < (-1)*alpha-h_posnega/sigma_stat
% % % sigma_stat^2*(z+alpha)^2
% func1_1 =   exp(-z^2/2)*(z-alpha)^2 ;
% a1 = int(func1_1, z, alpha-h_posnega,inf);
% func1_2 =  (z^2)*exp(-z^2/2);
% a2 = int(func1_2, z, (-1)*alpha-h_posnega, alpha-h_posnega);
% func1_3 = exp(-z^2/2)*(z+alpha)^2;
% a3 = int(func1_3, z, -inf, (-1)*alpha-h_posnega);
% % 
% % % case 0:  x = 0    0.872
% % % 1 : z>alpha
% % % sigma_stat^2*(z-alpha)^2
% % % 0: 0
% % 
% % % -1 : z< (-1)*alpha
% % % sigma_stat^2*(z+alpha)^2
% func2_1 =  exp(-z^2/2)*((z-alpha)^2) ;
% b1 = int(func2_1, z, alpha, inf);
% func2_2 = exp(-z^2/2)*(z^2);
% b2 = int(func2_2,z,(-1)*alpha,alpha);
% func2_3 = exp(-z^2/2)*(z+alpha)^2;
% b3 = int(func2_3, z, -inf, (-1)*alpha);
% % 
% % % case -1: x = -1    0.064
% % % 1: z>alpha+1/sigma_stat
% % % sigma_stat^2*(z-alpha)^2
% % 
% % % 0: (-1)*alpha+1/sigma_stat < z < alpha+1/sigma_stat
% % % 1
% % % -1: z < (-1)*alpha+1/sigma_stat
% % % % sigma_stat^2*(z+alpha)^2
% % 
% func3_1 = exp(-z^2/2)*(z-alpha)^2 ;
% c1 = int(func3_1, z, alpha+h_posnega, inf);
% func3_2 =  (z^2)*exp(-z^2/2);
% c2 = int(func3_2, z, (-1)*alpha+h_posnega, alpha+h_posnega);
% func3_3 = exp(-z^2/2)*(z+alpha)^2;
% c3 = int(func3_3, z, -inf, (-1)*alpha+h_posnega);
% 
% fMSE = ((a1+a2+a3)*epsilon/2+(b1+b2+b3)*(1-epsilon)+(c1+c2+c3)*epsilon/2)*(1/sqrt(2*pi));
% fMSE = epision*((sigma_e^2+sigma^2)*sigma_stat^2-sigma^2*sigma_e^2)/(sigma_e^2-sigma_stat^2);

% detection rate
% 1
% func =  1/sqrt(2*pi)*exp(-z^2/2);
% a1_dr  = int(func, z, (-1)*alpha-h_posnega/sigma_stat, alpha-h_posnega/sigma_stat);
% a1_conj = 1-vpa(a1_dr);
% 
% % 0
% a2_dr = int(func, z, (-1)*alpha, alpha);
% a2_conj = 1-vpa(a2_dr);
% 
% % -1
% a3_dr = int(func, z, (-1)*alpha+h_posnega/sigma_stat, alpha+h_posnega/sigma_stat);
% a3_conj = 1-vpa(a3_dr);
% 
% % 
% equ_dr = epsilon*a1_conj+(1-epsilon)*a2_conj;

syms z ;
syms h_posnega
% pi = 3.1416;
% lamda = 1e-2;
% tau = 3;
% case 1:  x = h_posnega    epision/2
% 1 : z>alpha-h_posnega/sigma_stat         % integral region
%  sigma_stat^2*(z-alpha)^2    % the function
% 0 : (-1)*alpha-h_posnega/sigma_stat < z < alpha-h_posnega/sigma_stat
% 1
% -1 : z < (-1)*alpha-h_posnega/sigma_stat
% % sigma_stat^2*(z+alpha)^2
func1_1 =   exp(-z^2/2)*(z-alpha)^2 ;
a1 = int(func1_1, z, alpha-h_posnega,inf);
func1_2 =  exp(-z^2/2)*(h_posnega^2);
a2 = int(func1_2, z, (-1)*alpha-h_posnega, alpha-h_posnega);
func1_3 = exp(-z^2/2)*(z+alpha)^2;
a3 = int(func1_3, z, -inf, (-1)*alpha-h_posnega);
% 
% % case 0:  x = 0    0.872
% % 1 : z>alpha
% % sigma_stat^2*(z-alpha)^2
% % 0: 0
% 
% % -1 : z< (-1)*alpha
% % sigma_stat^2*(z+alpha)^2
func2_1 =  exp(-z^2/2)*((z-alpha)^2) ;
b1 = int(func2_1, z, alpha, inf);
func2_2 = exp(-z^2/2)*0;
b2 = int(func2_2,z,(-1)*alpha,alpha);
func2_3 = exp(-z^2/2)*(z+alpha)^2;
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
func3_1 = exp(-z^2/2)*(z-alpha)^2 ;
c1 = int(func3_1, z, alpha+h_posnega, inf);
func3_2 =  (h_posnega^2)*exp(-z^2/2);
c2 = int(func3_2, z, (-1)*alpha+h_posnega, alpha+h_posnega);
func3_3 = exp(-z^2/2)*(z+alpha)^2;
c3 = int(func3_3, z, -inf, (-1)*alpha+h_posnega);

fMSE = ((a1+a2+a3)*epsilon/2+(b1+b2+b3)*(1-epsilon)+(c1+c2+c3)*epsilon/2)*(1/sqrt(2*pi));
eq = fMSE-mse_object;
%options = optimset('Display','iter');
h_scalar = solve(eq,h_posnega);
fMSE_num = subs(fMSE,'h_posnega',h_scalar);
% EqDR = subs(equ_dr, 'h_posnega', abs(double(h_scalar)));

% eq1 = sigma_stat^2-sigma_e^2*(sigma^2+fMSE/epision)/(sigma_e^2+sigma^2+fMSE/epision);
% % eq1 = lamda/tau-(sigma_e^2+sigma^2+fMSE/epision)/(sigma^2+fMSE/epision)*alpha*sigma_stat;
% eq2 = lamda-(sigma^2+fMSE/epision+sigma_e^2)/sigma_e^2*alpha*sigma_stat*(1-(sigma_e^2/(sigma_e^2+sigma^2+fMSE/epision))*equ_dr/epision);
% S = solve(eq1, eq2, alpha, sigma_stat);
% alpha_tru = double(S.alpha)
% sigma_stat_tru = double(S.sigma_stat)
% fMSE_tru = double(subs(fMSE,{alpha,sigma_stat},{double(S.alpha),double(S.sigma_stat)}))
% tau = lamda*(sigma^2+fMSE_tru/epision)/(alpha_tru*sigma_stat_tru*(sigma_e^2+sigma^2+fMSE_tru/epision))
 %[mse_empirical]= MSE_EMP(lamda,tau,noise_ratio)