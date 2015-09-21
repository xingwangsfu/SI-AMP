% This file is to test the application of GENP-AMP in Hybrid Multi-View
% Imaging System

clear all;
% close all;
row = 1024;
col = 768;
load Balloons_orig.mat;
load Balloons_interp.mat;
% img_origin=rgb2gray(img_origin); % rgb to grayscale
% img_interp=rgb2gray(img_interp);
[M,N] = size(img_origin);
sigma_trans = 1e3;  % variance of Gaussian noise added in the interpolated middle image
img_interp_noisy = double(img_interp)+ sqrt(sigma_trans)*randn(size(img_interp));% rgb to grayscale
x_mat = double(img_origin);
x_interp_mat = double(img_interp_noisy); %  interpolated middle image with Gaussian noise

% patch size 48 times 48 %
M_size = 48;
N_size = 48;
n_size = M_size*N_size;

dct_row = dct(eye(M_size));
dct_col = dct(eye(N_size));
Phi = kron(dct_col, dct_row); % DCT basis

delta_cand = [1/5];
sigma_SI_est_vector = cell(size(delta_cand));
sigma = 1e2; % variance of measurement noise added in the CS measurements
delta = delta_cand;
m_size = floor(delta*n_size);

% normalize the Gaussian measurement matrix %
Psi = randn(m_size,n_size);
A = Psi;
normalized_vector = sqrt(sum(A.*A,1));
normalized_matrix = repmat(normalized_vector,m_size,1);
A = A./normalized_matrix;
At = A';

% overlapping size 6 times 6 %
blocksize = [M_size N_size];
stepsize = [42 42];

if ~isa(A, 'function_handle')
    At1 = @(x) A'*x;
    A1 = @(x) A*x;
end
B = A*Phi';
Bt = B';
if ~isa(B, 'function_handle')
    Bt1 = @(x) B'*x;
    B1 = @(x) B*x;
end


%%%%% decompose a full size image into overlapping blocks %%%%%%%

% lastids contains the indices of the last block in each dimension
lastids = stepsize .* floor((size(x_mat)-blocksize)./stepsize) + 1;
dim_overlapp= floor((size(x_mat)-blocksize)./stepsize) + 1;
blocknum = prod(dim_overlapp);
ids_mtx = cell((dim_overlapp));
ids_mtx{1,1} = [1 1];
for ids_idx = 1:size(ids_mtx,1)
    for ids_idy = 1:size(ids_mtx,2)
        ids_mtx{ids_idx, ids_idy} = [(ids_idx-1)*stepsize(1)+1 (ids_idy-1)*stepsize(2)+1];
    end
end
ids_vec = reshape(ids_mtx, blocknum,1);

% matlabpool open 4 ; % open 4 pools to accelarate the computation
for i = 1 : blocknum
    ids = cell(2,1);
    ids{1} = [ids_vec{i}(1):ids_vec{i}(1)+M_size-1];
    ids{2} = [ids_vec{i}(2):ids_vec{i}(2)+N_size-1];
    block = x_mat(ids{1},ids{2}); % patchcs of original middle image
    theta = Phi*block(:);
    theta = theta(:);
    block = block(:);
    block_interp = x_interp_mat(ids{1},ids{2});  % patches of interpolated middle image
    theta_interp = Phi*block_interp(:);
    theta_interp = theta_interp(:);
    block_interp = block_interp(:);
    % get the CS measurements
    y = A1(block);
    y = y + sqrt(sigma)*randn(size(y));
    sigma_SI = norm(theta-theta_interp)^2/prod(size(theta));
   % [theta_ell1_ell1(:,i)] = CS_ell1_ell1(y, B, theta_interp, sigma);
    % Expectation maximization Gaussian Mixture model (EMGMAMP) %
%     [xhat, EMfin] = EMGMAMP(y, B);
%     [Xhat_v2, EMfin_v2] = EMGMAMP_v2(y, B,theta_interp); % generalized elastic net prior based EMGMAMP (GENP-EM)
%     theta_EM(:,i) = xhat;
%     theta_en_EM(:,i) = Xhat_v2;
%     
%     %%%% reconstruction algorithms %%%%
    [tmp_AMP, tmp_P_en_AMP, tmp_en_AMP, tmp_AMP_residual, sigma_SI_est, r_est] = patch_reconstruction (y, B, Bt, theta_interp, sigma_SI, 'Auto', theta, A1, At1);
    theta_AMP(:,i) =  ((tmp_AMP)); % AMP result
    theta_P_en_AMP(:,i) = ((tmp_P_en_AMP)); % Parameterless GENP-AMP
    theta_en_AMP(:,i) = ((tmp_en_AMP)); % GENP-AMP with known variance
    theta_AMP_residual(:,i) =  ((tmp_AMP_residual)); % AMP residual
%     
%     %%% denoise blocks by blocks
%     [tmp_denoise] = auto_soft_threshold(theta_interp,sigma_SI);
%     theta_denoise(:,i) = tmp_denoise; % soft-threshold denoising
%     percentage = 0.01;
%     [theta_ModiCS(:,i)] = ModiCS(y, B, theta_interp, percentage, sigma); % Modified CS
  
end


%%% reunite the patches into full size image %%%

Params.M = size(img_origin,1);
Params.N = size(img_origin,2);
Params.m = M_size;
Params.n = N_size;
Params.blocknum = blocknum;
Params.ids_vec = ids_vec;
Params.blocksize = blocksize;
Params.stepsize = stepsize;
M_dim = lastids(1)+blocksize(1)-1;
N_dim = lastids(2)+blocksize(2)-1;
Params.M_dim = M_dim;
Params.N_dim = N_dim;
img_origin_ref = img_origin(1:M_dim,1:N_dim);
[img_AMP,img_AMP_residual,img_en_AMP,img_P_en_AMP,img_EM,img_en_EM,img_denoise,img_ModiCS] = patch2image(theta_AMP,theta_AMP_residual,theta_en_AMP,theta_P_en_AMP,theta_EM,theta_en_EM, theta_denoise, theta_ModiCS, blocknum,Params);
PSNR_AMP = 20*log10(255/sqrt(mean((img_AMP(:)-img_origin_ref(:)).^2)));
PSNR_AMP_residual = 20*log10(255/sqrt(mean((img_AMP_residual(:)-img_origin_ref(:)).^2)));
PSNR_en_AMP = 20*log10(255/sqrt(mean((img_en_AMP(:)-img_origin_ref(:)).^2)));
PSNR_P_en_AMP = 20*log10(255/sqrt(mean((img_P_en_AMP(:)-img_origin_ref(:)).^2)));
PSNR_EM = 20*log10(255/sqrt(mean((img_EM(:)-img_origin_ref(:)).^2)));
PSNR_en_EM = 20*log10(255/sqrt(mean((img_en_EM(:)-img_origin_ref(:)).^2)));
PSNR_denoise = 20*log10(255/sqrt(mean((img_denoise(:)-img_origin_ref(:)).^2)));
PSNR_ModiCS = 20*log10(255/sqrt(mean((img_ModiCS(:)-img_origin_ref(:)).^2)));
matlabpool close;




