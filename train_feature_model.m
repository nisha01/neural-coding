function fea_mo = train_feature_model(im, var)
%TRAIN_FEATURE_MODEL
%   TODO Summary, detailed explanation and comment
%   This function was created by Umut Güçlü and last modified on 2 May 2014
%   at 21:15
%   Please visit http://www.ccnlab.net/ for more information

[im_size(1), im_size(2), im_num] = size(im);

for i = var.im_patch_num : -1 : 1
    
    im_patch(:, :, i) = im((1 : var.im_patch_size(1)) + randi(im_size(1) - (var.im_patch_size(1) + var.im_patch_off(1))) + var.im_patch_off(1), (1 : var.im_patch_size(2)) + randi(im_size(2) - (var.im_patch_size(2) + var.im_patch_off(2))) + var.im_patch_off(2), randi(im_num));
    
end

[fea_mo.neigh_func, fea_mo.fea_de] = tica(reshape(im_patch, var.im_patch_size(1) * var.im_patch_size(2), var.im_patch_num), var);

fea_mo.fea_de_size = var.im_patch_size;

end

function [neigh_func, fea_de] = tica(im_patch, var)

fea_num = var.lay_size(1) * var.lay_size(1);

fea_de = randn(fea_num, fea_num);
fea_de = (fea_de * fea_de') ^ -0.5 * fea_de;
fea_de = {fea_de};

min_neg_log_like_in = 1;

im_patch = bsxfun(@minus, im_patch, mean(im_patch, 2));

[eigenvec, eigenval] = svd(im_patch * im_patch' / var.im_patch_num);

eigenval = diag(eigenval);

whit_ma = diag(eigenval(1 : fea_num) .^ -0.5) * eigenvec(:, 1 : fea_num)';

im_patch = whit_ma * im_patch;

neigh_func = zeros(fea_num, var.lay_size(1), var.lay_size(2));

i = fea_num;

for j = var.lay_size(1) : -1 : 1
    
    for k = var.lay_size(1) : -1 : 1
        
        neigh_func(i, 1 : var.neigh_size(1), 1 : var.neigh_size(2)) = 1;
        
        neigh_func = circshift(neigh_func, [0 0 -1]);
        
        i = i - 1;
        
    end
    
    neigh_func = circshift(neigh_func, [0 -1, 0]);
    
end

neigh_func = reshape(neigh_func, fea_num, fea_num);

neigh_func_trans = neigh_func';

step_size = 1;

lay_size_times_im_patch_size = var.lay_size .* var.im_patch_size;

tic;

for i = var.it_num : -1 : 1
    
    fea_de_times_im_patch = fea_de{min_neg_log_like_in} * im_patch;
    
    gra = neigh_func_trans * (1 + neigh_func * fea_de_times_im_patch .^ 2) .^ -1 .* (2 * fea_de_times_im_patch) * im_patch';
    gra = gra - fea_de{min_neg_log_like_in} * gra' * fea_de{min_neg_log_like_in};
    
    step_size = arrayfun(@(co) co * step_size(min_neg_log_like_in), [0.5 1 2]);
        
    fea_de = cellfun(@(fea_de) (fea_de * fea_de') ^ -0.5 * fea_de, arrayfun(@(step_size) fea_de{min_neg_log_like_in} - step_size * gra, step_size, 'UniformOutput', false), 'UniformOutput', false);
    
    [min_neg_log_like, min_neg_log_like_in] = min(cellfun(@(fea_de) sum(sum(log(1 + neigh_func * (fea_de * im_patch) .^ 2))), fea_de));
    
    if toc >= 1
        
        disp(['iteration: ' num2str(i) ', negative log likelihood: ' num2str(round(min_neg_log_like)) ', step size: ' num2str(step_size(min_neg_log_like_in))]);
        
        fea = pinv(fea_de{min_neg_log_like_in} * whit_ma);
        fea = bsxfun(@minus, fea, min(fea));
        fea = bsxfun(@rdivide, fea, max(fea));
        
        imagesc(col2im(fea, var.im_patch_size, lay_size_times_im_patch_size, 'distinct'));
        
        colormap(gray);
        
        axis off;
        
        title('feature detectors');
        
        drawnow;
        
        tic;
        
    end
    
end

fea_de = fea_de{min_neg_log_like_in} * whit_ma;

end