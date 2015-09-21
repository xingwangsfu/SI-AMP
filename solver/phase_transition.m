clear all;

syms alpha;
psi = exp(-1*alpha^2/2)/sqrt(2*pi);

Phi = int(psi, -inf, -alpha);

var_epision = 1;

rho_temp = (1 - alpha*Phi/psi)*1;

LR_factor = 4;
delta = 0.4;
delta_LR = delta;

for i = 1:length(delta)
   
%     eqs = (2*psi/((alpha+2*(psi-alpha*Phi))))-LR_factor*delta(i);
%     alpha_num = solve(eqs,alpha);
%     rho(i) = double(subs(rho_temp,'alpha', abs(double(alpha_num))));
    
    eqs_LR = 2*psi/(alpha+2*(psi-alpha*Phi))*1-delta_LR(i);
    alpha_num_LR = solve(eqs_LR,alpha);
    rho_LR(i) = double(subs(rho_temp,'alpha',abs(double(alpha_num_LR))));
    [alpha,lambda,rho_DMM] = optimal_threshold_DMM(delta(i));
    rho_DMM_LR(i) = double(rho_DMM);
    
    i
end


