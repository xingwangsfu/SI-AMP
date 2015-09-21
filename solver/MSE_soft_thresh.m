function mse_soft_thresh = MSE_soft_thresh(alpha, sigma_stat, h_posnega, epsilon)
% syms sigma_stat;

syms h_posnega;
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
% % sigma_stat^2*(z+alpha)^2;
func1_1 =   1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z-alpha)^2 ;
a1 = int(func1_1, z, alpha-h_posnega/sigma_stat,inf);
func1_2 =  1/sqrt(2*pi)*exp(-z^2/2)*h_posnega^2;
a2 = int(func1_2, z, (-1)*alpha-h_posnega/sigma_stat, alpha-h_posnega/sigma_stat);
func1_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
a3 = int(func1_3, z, -inf, (-1)*alpha-h_posnega/sigma_stat);
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
c1 = int(func3_1, z, alpha+h_posnega/sigma_stat, inf);
func3_2 =  1/sqrt(2*pi)*exp(-z^2/2)*h_posnega^2;
c2 = int(func3_2, z, (-1)*alpha+h_posnega/sigma_stat, alpha+h_posnega/sigma_stat);
func3_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
c3 = int(func3_3, z, -inf, (-1)*alpha+h_posnega/sigma_stat);

mse_soft_thresh = (a1+a2+a3)*(epsilon/2)+(b1+b2+b3)*(1-epsilon)+(c1+c2+c3)*(epsilon/2);