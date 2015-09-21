function [theta_ModiCS_vec] =ModiCS(y, B, theta_interp, percentage, sigma)
    n_size = prod(size(theta_interp));
    m_size = prod(size(y));
 %   
    full_idx = [1:n_size];
    support_idx = support_detection(theta_interp, percentage);
    support_idx_complementary =  setdiff(full_idx, support_idx);
    %  para_cvx = delta_cand*(percentage*sum(theta_interp.^2) + (1-percentage)*sigma_SI) + sigma;
     %para_cvx = sigma;
    
    % incorporate cvx command into a function
    cvx_begin quiet
        variable theta_tmp(n_size);
        minimize (norm(theta_tmp(support_idx_complementary),1))
    subject to
      %  y == B*theta_ModiCS
        norm(y-B*theta_tmp(:),2) <= sqrt(m_size*sigma)
    cvx_end
    
    theta_ModiCS_vec = theta_tmp;
