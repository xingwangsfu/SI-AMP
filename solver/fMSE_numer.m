function fMSE = fMSE_numer(h_posnega,alpha,epsilon)
syms z ;
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