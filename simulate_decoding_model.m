function im_hat = simulate_decoding_model(fea_mo, im, vox, vox_mo, vox_num)
%SIMULATE_DECODING_MODEL
%   TODO Summary, detailed explanation and comment
%   This function was created by Umut Güçlü and last modified on 2 May 2014
%   at 21:15
%   Please visit http://www.ccnlab.net/ for more information

im_num = size(im, 3);

vox_hat = simulate_voxel_model(simulate_feature_model(fea_mo, im), vox_mo);

vox_trans     = vox';
vox_hat_trans = vox_hat';

for i = im_num : -1 : 1
    
    disp(['image: ' num2str(i)]);
    
    im_in = i ~= 1 : im_num;
    
    [~, de_rho_in] = sort(diag(corr(vox(im_in, :), vox_hat(im_in, :))), 'descend');
    
    vox_in = de_rho_in(1 : vox_num);
    
    [~, im_hat(i)] = max(corr(vox_trans(vox_in, i), vox_hat_trans(vox_in, :)));
    
end

end

