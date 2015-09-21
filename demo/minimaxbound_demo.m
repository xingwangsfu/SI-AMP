% this document is to check the minimax risk of SI-LASSO, corresponding to
% the first experiment in the papaer
clear all;
varsigma_e_square_cand = 4; % the ratio between variance of prediction 
                             %error and variance of measurement noise
delta = 0.1;
c = 0;
d = 0.02; % the nearly least favorable distribution for scalar 
          % soft-thresholding denoiser
sigma_w = 1; % variance of measurement noise
[epsilon] =1.9*delta;
N = 2000;
sigma_SI_est_vector= cell(size(varsigma_e_square_cand),1);
sigma_SI_est_vector_fake = cell(size(varsigma_e_square_cand),1);
epsilon_i = epsilon;   % the sparsity rate

for var_index = 1:length(varsigma_e_square_cand)
    varsigma_e_square = varsigma_e_square_cand(var_index);
    sigma_SI = varsigma_e_square*sigma_w;  % the variance of prediction error 
                                           % (generalized elastic net prior)
    
    for delta_index = 1:length(delta)
        sigma = sigma_w;
        i = delta_index;
        syms alpha;
        
        % get the minimax MSE for scalar case
        [alpha,alpha_num] = optimal_threshold(epsilon_i); % to get the minimax
                                                          % optimal thresholding parameter in soft-thresholding function
        psi = exp(-1*(alpha^2) /2)/sqrt(2*pi);
        Phi = int(psi, -inf, -1*alpha);
        mse_scalar_num = epsilon_i*(1+alpha^2)+(1-epsilon_i)*(2*(1+alpha^2)*Phi-2*alpha*psi);
        mse_scalar_simp = 2*psi/(alpha+2*psi-2*alpha*Phi);
        mse_scalar_simp = double(subs(mse_scalar_simp,'alpha',((alpha_num))));
        mse_scalar = double(subs(mse_scalar_num, 'alpha', ((alpha_num))));
        alpha_vector(i) = alpha_num;
        M_scalar(i) =  mse_scalar;
         
        % the minimax MSE for AMP
        %         m_scalar_object = double((1-d)*mse_scalar)
        %         MSE_classical(i) =  m_scalar_object/(1- m_scalar_object/delta(delta_index));
        %
        %         NPI(i) = 1 + MSE_classical(i)/delta(delta_index);
        %        
        
        % the nearly least-favorable 3-point distribution for minimax MSE
        [h_scalar] = solve_h_posnega(epsilon_i, (1-d)*M_scalar(i), alpha_vector(i));
        h_posnega(i) = double(h_scalar);
        
        
        % the minimax MSE for GENP-AMP ( MSE_SI_vector ) and its correponding nearly
        % least-favorable 3-point distribution ( h_SI )
        G_temp_exact = sqrt((delta(delta_index)*varsigma_e_square + delta(delta_index) -(1-c)*varsigma_e_square*M_scalar(i))^2+...
            4*(1-c)*delta(delta_index)*varsigma_e_square*M_scalar(i))+(delta(delta_index)*varsigma_e_square+delta(delta_index)-(1-c)*varsigma_e_square*M_scalar(i));
        G_temp_appro = sqrt((delta(delta_index)*varsigma_e_square + delta(delta_index) -(1-d)*varsigma_e_square*M_scalar(i))^2+...
            4*(1-d)*delta(delta_index)*varsigma_e_square*M_scalar(i))+(delta(delta_index)*varsigma_e_square+delta(delta_index)-(1-d)*varsigma_e_square*M_scalar(i));
        MSE_SI_exact = M_scalar(i)*2*(1-c)*delta(delta_index)*varsigma_e_square/G_temp_exact;
        MSE_SI_appro = M_scalar(i)*2*(1-d)*delta(delta_index)*varsigma_e_square/G_temp_appro;
        MSE_SI_vector(i) = double(MSE_SI_appro);
        NPI_SI = sigma*(varsigma_e_square*(1+MSE_SI_exact/delta(delta_index))/(varsigma_e_square+1+MSE_SI_exact/delta(delta_index)));
        fNPI_SI = sigma*(varsigma_e_square*(1+MSE_SI_appro/delta(delta_index))/(varsigma_e_square+1+MSE_SI_appro/delta(delta_index)));
        h_SI(i) = h_posnega(i)*sqrt(NPI_SI);
        mu_stat = (1+MSE_SI_appro/delta(delta_index))/varsigma_e_square;
        EqDR_SI = (Equ_DR(sqrt(fNPI_SI), alpha_vector(i), h_SI(i), epsilon_i));
        
        % the corrsponding parameters lambda and tau for  LASSO like formulate
        lambda_minimax_SI(i) = (1+mu_stat)*alpha_vector(i)*sqrt(fNPI_SI)*(1-1/(1+mu_stat)*EqDR_SI/delta(delta_index));
        tau_minimax_SI(i) = mu_stat*(1-1/(1+mu_stat)*EqDR_SI/delta(delta_index));
        
        
        for index = 1:20
            
            % generate the sparse signal following 3-point distribution
            x = zeros(N,1);
            rand_state = rand(N,1);
            x(rand_state<=epsilon_i/2) = h_SI(i);
            x(rand_state>=1-epsilon_i/2) = -h_SI(i);
            C_x = cov(x);
            Mu_x = mean(x);
            x_SI = x + sqrt(sigma_SI)*randn(size(x));
            
            % generat the Gaussian measurement matrix
            M = floor(delta(delta_index)*N);
            A = randn(M,N);
            normalized_vector = sqrt(sum(A.*A,1));
            normalized_matrix = repmat(normalized_vector,M,1);
            A = A./normalized_matrix;
            AT = A';
            Y = A*x + sqrt(sigma)*randn(M,1);
            [ err1 ] = mmwrite( 'A_matrix.mtx',A);
            % denoising the prediction (GENP) with soft-thresholding
            % directly
            x_denoise = wthresh(x_SI, 's', double(alpha_num)*sqrt(sigma_SI));
            eMSE_DN(index) = norm(x_denoise - x)^2/length(x);
            
            
            % GENP-AMP, minimax thresholding rule
            Params.x_SI = x_SI; %  side information
            Params.T = 60; % number of iterations
            Params.tol = 1e-8; % tolerance
            Params.sigma_SI = sigma_SI; % variance of prediction error
            Params.mode = 'Minimax'; % two modes here: auto and minimax, auto for parameterless AMP, minimax for DMM
            Params.SIorNot = 'Yes'; % two modes: SI exists or not
            [alpha_rec] = SI_AMP(Y, A, AT, N,double(alpha_num), Params );
            eMSE(delta_index,index) = norm(alpha_rec-x)^2/N;
            
            
            % LMMSE (linear minimum mean squared error): assuming x follows
            % Gaussian distribution
            A_extend = [A;eye(N)];
            Y_extend = [Y;x_SI];
            C_w = sigma_w*eye(M,M);
            C_e = sigma_SI*eye(N,N);
            C_w = [C_w;zeros(N,M)];
            C_e = [zeros(M,N);C_e];
            C_z = [C_w C_e];
            x_LMMSE = C_x*A_extend'*((A_extend*C_x*A_extend'+C_z)\(Y_extend-A_extend*(Mu_x*ones(N,1))))+Mu_x;
            eMSE_LS(delta_index,index) = norm(x_LMMSE-x)^2/N;
            %
            
            %  save the date into mtx form, served as the input to OWLQN
            %  solver
            %             [ err1 ] = mmwrite( 'A_matrix.mtx',A); % %
            %             [ err2] = mmwrite( 'y_matrix.mtx',Y); %
            %             [ err3 ] = mmwrite( 'SI_matrix.mtx',x_SI); %
            %             [err4 ] = mmwrite('x_matrix.mtx',x);
         
        end
    end
end

