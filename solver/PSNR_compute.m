function [PSNR_AMP, PSNR_para_en_AMP, PSNR_en_AMP, PSNR_AMP_residual] = PSNR_compute(theta_AMP, theta_parameterless_en_AMP, theta_en_AMP, theta_AMP_residual, theta)
n_size = prod(size(theta));

PSNR_AMP = 20*log10(255/(sqrt(norm(theta_AMP-theta)^2/n_size)));

PSNR_para_en_AMP = 20*log10(255/sqrt(norm(theta_parameterless_en_AMP-theta)^2/n_size));

PSNR_en_AMP = 20*log10(255/sqrt(norm(theta_en_AMP-theta)^2/n_size));

PSNR_AMP_residual = 20*log10(255/sqrt(norm(theta_AMP_residual-theta)^2/n_size));