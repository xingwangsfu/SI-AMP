function [alpha,alpha_num] = optimal_threshold(epsilon)   

% this function is to compute the optiam threhsolding parameter given the
% sparsity epsilon
   syms alpha;
    psi = exp(-1*(alpha^2) /2)/sqrt(2*pi);
    
    Phi = int(psi, -inf, -1*alpha);
    
    eqs = (2*psi-2*alpha*Phi)/(alpha+2*psi-2*alpha*Phi)-epsilon;
    % eqs = alpha*epsilon+(1-epsilon)*(2*alpha*Phi-2*psi);
    
    alpha_num = solve(eqs,alpha);
    
    alpha_num = abs(double(alpha_num));