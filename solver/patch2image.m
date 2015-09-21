function [img_AMP,img_AMP_residual,img_en_AMP,img_P_en_AMP,img_EM,img_en_EM,img_denoise,img_ModiCS] = patch2image(theta_AMP,theta_AMP_residual,theta_en_AMP,theta_P_en_AMP,theta_EM,theta_en_EM, theta_denoise, theta_ModiCS, blocknum,Params)

blocksize = Params.blocksize;
block_row = blocksize(1);
block_col = blocksize(2);

for i = 1:blocknum
    tmp_AMP = idct2(reshape(theta_AMP(:,i),block_row,block_col));
    tmp_AMP_residual =  idct2(reshape(theta_AMP_residual(:,i),block_row,block_col));
    tmp_en_AMP =  idct2(reshape(theta_en_AMP(:,i),block_row,block_col));
    tmp_P_en_AMP =  idct2(reshape(theta_P_en_AMP(:,i),block_row,block_col));
    tmp_EM =  idct2(reshape(theta_EM(:,i),block_row,block_col));
    tmp_EM_v2 =  idct2(reshape(theta_en_EM(:,i),block_row,block_col));
    tmp_denoise =  idct2(reshape(theta_denoise(:,i),block_row,block_col));
    tmp_ModiCS =  idct2(reshape(theta_ModiCS(:,i),block_row,block_col));
    x_AMP(:,i) = tmp_AMP(:);
    x_AMP_residual(:,i) = tmp_AMP_residual(:);
    x_en_AMP(:,i) = tmp_en_AMP(:);
    x_P_en_AMP(:,i) = tmp_P_en_AMP(:);
    x_EM(:,i) = tmp_EM(:);
    x_en_EM(:,i) = tmp_EM_v2(:);
    x_denoise(:,i) = tmp_denoise(:);
    x_ModiCS(:,i) = tmp_ModiCS(:);
end

[img_AMP] = patch_unite(x_AMP,Params);
[img_en_AMP] = patch_unite(x_en_AMP,Params);
[img_P_en_AMP] = patch_unite(x_P_en_AMP,Params);
[img_AMP_residual] = patch_unite(x_AMP_residual,Params);
[img_EM] = patch_unite(x_EM,Params);
[img_en_EM] = patch_unite(x_en_EM,Params);
[img_denoise] = patch_unite(x_denoise,Params);
[img_ModiCS] = patch_unite(x_ModiCS,Params);