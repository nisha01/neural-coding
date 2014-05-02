function im_fea = simulate_feature_model(fea_mo, im)
%SIMULATE_FEATURE_MODEL
%   TODO Summary, detailed explanation and comment
%   This function was created by Umut Güçlü and last modified on 2 May 2014
%   at 19:42
%   Please visit http://www.ccnlab.net/ for more information

im_num = size(im, 3);

for i = im_num : -1 : 1
    
    disp(['image: ' num2str(i)]);
    
    im_fea(:, :, i) = log(1 + fea_mo.neigh_func * (fea_mo.fea_de * im2col(im(:, :, i), fea_mo.fea_de_size, 'distinct')) .^ 2);
    
end

end