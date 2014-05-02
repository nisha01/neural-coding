function lamb_hat = solve_df_for_lambda(df, s, var)
%SOLVE_DF_FOR_LAMBDA
%   TODO Summary, detailed explanation and comment
%   This function was created by Umut Güçlü and last modified on 2 May 2014
%   at 19:42
%   Please visit http://www.ccnlab.net/ for more information

df = sort(df, 'ascend');

guess   = @(df, lamb, N, s) max(1 / mean(s .^ -2) * (N - df) / df, lamb);
f       = @(df, lamb, s)  sum(s .^ 2 ./ (s .^ 2 + lamb)) - df;
f_prime = @(lamb, s) -sum(s .^ 2 ./ (s .^ 2 + lamb) .^ 2);
iterate = @(df, lamb, s) max(0, lamb - f(df, lamb, s) / f_prime(lamb, s));

for i = var.df_num : -1 : 1
    
    disp(['df: ' num2str(i)]);
    
    if i == var.df_num
        
        lamb_hat(i) = guess(df(i), 0 , length(s), s);
        
    else
        
        lamb_hat(i) = guess(df(i), lamb_hat(i + 1), length(s), s);
        
    end
    
    while abs(f(df(i), lamb_hat(i), s)) > var.df_tol
        
        lamb_hat(i) = iterate(df(i), lamb_hat(i), s);
        
    end
    
end

end