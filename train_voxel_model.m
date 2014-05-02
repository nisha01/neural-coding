function vox_mo = train_voxel_model(im_fea, var, vox)
%TRAIN_VOXEL_MODEL
%   TODO Summary, detailed explanation and comment
%   This function was created by Umut Güçlü and last modified on 2 May 2014
%   at 21:15
%   Please visit http://www.ccnlab.net/ for more information

[fea_num, block_num, im_num] = size(im_fea);

im_fea = reshape(im_fea, fea_num * block_num, im_num);

im_fea_mu = mean(im_fea, 2);

im_fea = bsxfun(@minus, im_fea, im_fea_mu);

[U, S, V] = svd(im_fea, 'econ');

s = diag(S);

s_squared = s .^ 2;

vox_mo.lamb_hat = solve_df_for_lambda(linspace(1, im_num, var.df_num), s, var);

V_trans_times_vox = V' * vox;

for i = var.df_num : -1 : 1
    
    disp(['lambda: ' num2str(i)]);
    
    s_square_plus_lamb = s_squared + vox_mo.lamb_hat(i);
    
    tem_Be_hat(:, :, i) = U * diag(s ./ s_square_plus_lamb) * V_trans_times_vox;
        
    GCV(:, i) = mean(((vox - im_fea' * tem_Be_hat(:, :, i)) / (1 - mean(s_squared ./ s_square_plus_lamb))) .^ 2);
    
end

Be_hat_num = size(tem_Be_hat, 2);

[~, min_GCV_in] = min(GCV, [], 2);

for i = Be_hat_num : -1 : 1
    
    disp(['Beta: ' num2str(i)]);
    
    vox_mo.Be_hat(:, i) = tem_Be_hat(:, i, min_GCV_in(i));
        
end

vox_mo.Be_hat(end + 1, :) = mean(vox) - im_fea_mu' * vox_mo.Be_hat;

end