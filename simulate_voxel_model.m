function vox_hat = simulate_voxel_model(im_fea, vox_mo)
%SIMULATE_VOXEL_MODEL
%   TODO Summary, detailed explanation and comment
%   This function was created by Umut Güçlü and last modified on 2 May 2014
%   at 19:42
%   Please visit http://www.ccnlab.net/ for more information

[fea_num, block_num, im_num] = size(im_fea);

im_fea = reshape(im_fea, fea_num * block_num, im_num);

im_fea(end + 1, :) = 1;

vox_hat = im_fea' * vox_mo.Be_hat;

end