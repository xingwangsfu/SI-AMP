function [theta_AMP, theta_parameterless_en_AMP, theta_en_AMP, theta_AMP_residual, sigma_SI_est, r_est ] = patch_reconstruction (y, A1, At1, theta_interp, sigma_SI, mode, theta, Psi1,Psit1t)


[m_size] = size(y,1);
n_size = size(theta_interp,1);
delta = m_size/n_size;

% Parameterless GENP-AMP %
for est_index = 1:2
    if est_index == 1
        [theta_AMP, r_tilde_real] =  SI_AMP(y, A1, At1, n_size, 0); % perform parameterless AMP
    else
        %   tic;
        Params.x_SI = theta_interp; % no side information
        Params.T = 30; % number of iterations
        Params.tol = 1e-8; % tolerance
        Params.sigma_SI = sigma_SI_est;
        Params.mode = 'Auto'; % two modes here: auto and minimax, auto for parameterless AMP, minimax for DMM
        Params.SIorNot = 'Yes'; % two modes: SI exists or not
        [theta_parameterless_en_AMP,  r_tilde_real_SI_en_AMP] =  SI_AMP(y, A1, At1, n_size, 0, Params);
    end
    sigma_SI_est = (norm(theta_interp - theta_AMP)^2 - (r_tilde_real(end)))/n_size;
    r_est = r_tilde_real(end)/n_size;
    if est_index == 1
        sigma_SI_est_copy = sigma_SI_est;
    end
end

% GENP-AMP with true variance of prediction error %
Params.x_SI = theta_interp;
Params.T = 60;
Params.tol = 1e-8;
Params.sigma_SI = sigma_SI;
Params.mode = 'Auto';
Params.SIorNot = 'Yes';
[theta_en_AMP] = SI_AMP(y, A1, At1, n_size, 0, Params);

% Reconstruting the predition error with AMP
y_error = y - A1*theta_interp;
theta_interp_without = sqrt(1e14)*randn(n_size,1);
sigma_interp_without = 1e14;
[theta_error, r_tilde_real] = SI_AMP(y_error, A1, At1, n_size, 0);
theta_AMP_residual = theta_error + theta_interp;



