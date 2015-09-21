function [theta] = CS_ell1_ell1(y, B, theta_interp, sigma)

n_size = prod(size(theta_interp));
m_size = prod(size(y));

cvx_begin quiet
   variable theta_tmp(n_size);
   minimize (norm(theta_tmp,1)+norm(theta_tmp-theta_interp,1))
subject to
    norm(y-B*theta_tmp(:),2) <= sqrt(m_size*sigma)
cvx_end


theta = theta_tmp;