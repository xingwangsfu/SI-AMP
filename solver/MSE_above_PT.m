% for mse(0,tau) < delta
clear all;
syms alpha_num;
syms z;
sigma_stat = 1;
delta = 0.2;
varsigma_e_square = 1;

% func2_1 =   1/sqrt(2*pi)*exp(-z^2/2)*(sigma_stat^2)*((z-alpha)^2) ;
% b1 = int(func2_1, z, alpha, inf);
% func2_2 = 0;
% b2 = 0;
% func2_3 = 1/sqrt(2*pi)*exp(-z^2/2)*sigma_stat^2*(z+alpha)^2;
% b3 = int(func2_3, z, -inf, (-1)*alpha);
%
% equ = b1+b2+b3-(1+1/varsigma_e_square)*delta;
%
% alpha_scalar = solve( equ, alpha);
%
epsilon = delta*0.28;

[beta,alpha] =  optimal_threshold_DMM(delta);

N = 2000;
gamma_cand = [ 0.99];
sigma = 4;
for i = 1 : length(gamma_cand)
    gamma = gamma_cand(i);
    % worst case MSE
    % syms alpha;
    psi = exp(-1*(alpha_num^2) /2)/sqrt(2*pi);
    
    Phi = int(psi, -inf, -1*alpha_num);
    
    psi_num = double(subs(psi,'alpha_num', alpha));
    
    Phi_num = double(subs(Phi, 'alpha_num', alpha));
    
    MSE_posnega = epsilon*(1+alpha^2)+(1-epsilon)*(2*(1+alpha^2)*Phi_num-2*alpha*psi_num);
    
    d = 1 - (gamma)*delta/MSE_posnega;
    
    m_scalar_object = (gamma)*delta;
    
    [h_scalar] = solve_h_posnega(epsilon, m_scalar_object, alpha);
    
    Q_gamma_varsigma = sqrt((1-gamma)^2+4*(1+gamma*varsigma_e_square)/varsigma_e_square^2)+(1-gamma);
    
    fMSE_SI = gamma*delta/(1-gamma);
    
    fNPI_SI = 1/(1-gamma);
    
    MSE(i) = fMSE_SI;
    
    NPI(i) = sigma*fNPI_SI;
    
    h_posnega(i) = abs(double(h_scalar));
    
    h_SI(i) = h_posnega(i)*sqrt(NPI(i));
    
    mu_stat = (1+fMSE_SI/delta)/varsigma_e_square;
    
    EqDR_SI = double(Equ_DR(sqrt(fNPI_SI), alpha, h_SI(i), epsilon));
    
    lamda_minimax_SI(i) = (1+mu_stat)*alpha*sqrt(fNPI_SI)*(1-1/(1+mu_stat)*EqDR_SI/delta);
    
    tau_minimax_SI(i) = mu_stat*(1-1/(1+mu_stat)*EqDR_SI/delta);
    
    psi = exp(-1*(alpha_num^2) /2)/sqrt(2*pi);
    Phi = int(psi, -inf, -1*alpha_num);
    [beta,lambda]  = optimal_threshold_DMM(2*delta);
    [beta,lambda_LR] = optimal_threshold(2*epsilon);
    epsilon = epsilon;
    mse_scalar_num = 2*epsilon*(1+alpha_num^2)+(1-2*epsilon)*(2*(1+alpha_num^2)*Phi-2*alpha_num*psi);
    mse_scalar_simp = 2*psi/(alpha_num+2*psi-2*alpha_num*Phi);
    mse_scalar = double(subs(mse_scalar_simp,'alpha_num',lambda_LR));
    fMSE_LR = mse_scalar/(1-mse_scalar/(2*delta));
    t_HR = zeros(100,1);
    t_LR = zeros(100,1);
    parfor index_x = 1:100
        x_ref = zeros(N/2,1);
        rand_state = rand(N/2,1);
        x_ref(rand_state<=2*epsilon/2) = h_SI(i);
        x_ref(rand_state>=1-2*epsilon/2) = -h_SI(i);
        e_appro = randn(N/2,1);
        % measument vector
        x_ref = [x_ref;zeros(N/2,1)];
        % x = [x;e_appro];
        M = floor(delta*N);
        A = randn(M,N);
        normalized_vector = sqrt(sum(A.*A,1));
        normalized_matrix = repmat(normalized_vector,M,1);
        A = A./normalized_matrix;
        y = A*x_ref;
        w = sqrt(sigma)*randn(M,1);
        y = y+w;
        AT = A';
        % side information
        tic;
        [alpha_rec] = SI_AMP(y, A, AT, N, alpha, x_ref, 1 );
        t_HR(index_x) = t_HR(index_x)+toc;
        eMSE(index_x) = norm(alpha_rec-x_ref)^2/N;
        
        A_LR = A(:,1:N/2);
        AT_LR = A_LR';
        tic;
        [alpha_rec_LR,xi] = SI_AMP(y, A_LR, AT_LR, N/2, lambda_LR, x_ref(1:N/2), 1 );
        t_LR(index_x) = t_LR(index_x)+toc;
        eMSE_LR(index_x) = norm(alpha_rec_LR-x_ref(1:N/2))^2/(N/2);
        
    end
end
