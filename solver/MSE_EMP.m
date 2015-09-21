% this document is to justify the prediction performance of
% side-information based AMP algorithm
% function [mse_empirical]= MSE_EMP(lamda,tau,noise_ratio)
% source
clear all;
N = 2000;
% mse_theoretical = 0;
% mse_epirical = 0;

delta = 0.25;
sigma = sqrt(0.2);
ratio_noise = 2;
sigma_e =ratio_noise*sqrt(0.2);
lamda_cand = 0.11:0.1:2;
for index_lamda = 1:length(lamda_cand)
    lamda = lamda_cand(index_lamda);
    [fMSE_tru] = get_lamda_tau(sigma, sigma_e, epsilon, lambda);
    %     MSE_theoretical(index_lamda) = fMSE_tru;
    %     tau_theoretical(index_lamda) = tau;
    for index_x = 1:20
        x = zeros(N,1);
        rand_state = rand(N,1);
        epsilon = 0.128;
        x(rand_state<=0.128) = randn(length(find(rand_state<=0.128)),1);
        
        
        % measument vector
        
        M = floor(delta*N);
        A = randn(M,N);
        normalized_vector = sqrt(sum(A.*A,1));
        normalized_matrix = repmat(normalized_vector,M,1);
        A = A./normalized_matrix;
        y = A*x;
        w = sigma*randn(M,1);
        y = y+w;
        % side information
        %   opt
        prediction_e = sigma_e*randn(N,1);
        x_SI = x + prediction_e;
        mse_SI(index_x) = norm(x_SI-x)^2/N;
        %
        %         [ err1 ] = mmwrite( 'A_matrix.mtx',A);
        %
        %         [ err2 ] = mmwrite( 'y_matrix.mtx',y);
        %         [ err3 ] = mmwrite( 'SI_matrix.mtx',x_SI);
        %         [err4 ] = mmwrite( 'x_matrix.mtx',x);
        % [x_rec, r_tilde_real,x_rec_unthresh, xi]= SI_AMP(y, A, A', M, N, x_SI, sigma^2, sigma_e^2, alpha, epsion, T, tol, delta,x,'auto','Not');
        [alpha_rec] = SI_AMP(Y, A, AT, N, alpha_num, x, 1 );
            eMSE(delta_index,index) = norm(alpha_rec-x)^2/N;
        % compute the stationary point
        
        % alpha_cand = 0.001;
        
        %
        % lamda = 1;
        % tau =  1.1605;
        % GPSR
        %
        %         cvx_begin quiet
        %         variable theta(N);
        %         minimize (1/2* sum((A*theta-y).^2)+lamda*norm(theta,1)+tau/2*sum((x_SI-theta).^2));
        %         cvx_end
        %         %          At = A';
        %         %         [theta,alp_debias,objective,times,debias_start,mses]= ...
        %         %         GPSR_BB_Modi(y,A,lamda,tau,x_SI,...
        %         %         'AT', At,'Debias',0,'Initialization',0,...
        %         %         'StopCriterion',1,'ToleranceA',0.0001,'ToleranceD',0.0001);
        %
        %         mse_empirical(index_lamda,index_x) = norm(theta-x)^2/N;
    end
end

% mse_theoretical = mse_theoretical/20
% mse_epirical = mse_epirical/20
