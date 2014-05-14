function G = Gabor_wavelet_pyramid(b, FOV, gamma, lambda, sigma, theta)
%GABOR_WAVELET_PYRAMID Summary of this function goes here
%   Detailed explanation goes here

length_of_FOV    = length(FOV);
length_of_theta  = length(theta);
length_of_x_0    = length_of_FOV ./ lambda;

G    = cell(2, 1);
G{1} = zeros(length_of_FOV, length_of_FOV, sum(length_of_theta * length_of_x_0 .^ 2));
G{2} = zeros(size(G{1}));

i = 1;

for j = 1 : length(lambda)
    
    x_0 = linspace(FOV(1), FOV(end), length_of_x_0(j));
    
    for k = 1 : length_of_theta
        
        for l = 1 : length_of_x_0(j)
            
            for m = 1 : length_of_x_0(j)
                
                G{1}(:, :, i) = Gabor_wavelet(b, FOV, gamma, lambda(j), 0     , sigma, theta(k), x_0(l), x_0(m));
                G{2}(:, :, i) = Gabor_wavelet(b, FOV, gamma, lambda(j), pi / 2, sigma, theta(k), x_0(l), x_0(m));
                
                i = i + 1;
                
            end
            
        end
        
    end
    
end


end

