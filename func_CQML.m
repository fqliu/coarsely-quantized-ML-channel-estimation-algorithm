function h = func_CQML(thres_l, thres_r, y_lowb, X, Nr)
%% Input
% thres_l : lower bound vector (thres_l = [l_R^T, l_I^T]^T, where l_R, l_I can
% be found in eq. (6) in [2])

% thres_r : upper bound vector (thres_r = [u_R^T, u_I^T]^T, where u_R, u_I can
% be found in eq. (6) in [2])

% y_lowb : low-bit quantized value (real-valued vector)
% X : pilot signal
% Nr : receive antenna number
%% Output
% h : estimated real-valued channel vector
%% Reference: [2] F. Liu, X. Shang, Y. Cheng and G. Zhang,¡°Computationally effcient maximum
%% likelihood channel estimation for coarsely quantized massive MIMO systems,¡±
%% IEEE Communications Letters, vol. 26, no. 2, Feb. 2022. 

[Nt, K] = size(X);

y_lowb = y_lowb(1:K*Nr) + 1i * y_lowb(K*Nr+1:2*K*Nr);
Y_lowb = reshape(y_lowb,Nr,K); 

t_l = thres_l(1:K*Nr) + 1i * thres_l(K*Nr+1:2*K*Nr);
T_l = reshape(t_l,Nr,K); 
T_lr = real(T_l);
T_li = imag(T_l);

t_r = thres_r(1:K*Nr) + 1i * thres_r(K*Nr+1:2*K*Nr);
T_r = reshape(t_r,Nr,K);
T_rr = real(T_r);
T_ri = imag(T_r);

% toc
Xt = X'/(X*X');  %% X^dag in eq. (20)


%% Initialization
H_hat_pre = Y_lowb * Xt;  %% initialize the channel matrix estimate
% Y_lowb - H_hat_pre * X
sigma_pre = norm(Y_lowb - H_hat_pre * X,'fro')/sqrt(2*K*Nr);
gamma_pre = 1 / sigma_pre;

if K == Nt
    sigma_pre = sqrt(mean(abs(y_lowb).^2)/2);
    gamma_pre = 1 / sigma_pre;
end


epsilon_1 = 1e-5;
epsilon_2 = 1e-4;

I_max_1 = 50; %% maximum iteration number in the outer cyclic iterations
I_max_2 = 10; %% maximum iteration number in the inner MM iterations

H_hat_pre = H_hat_pre * gamma_pre;
fx_cdf = @(x) 0.5 * erfc(-1/sqrt(2) * x );

relative_dif_1 = 1;
iter_num_1 = 0;
while(relative_dif_1 > epsilon_1) %% Outer cyclic iteration
%     tic
   iter_num_1 = iter_num_1 + 1;
   if iter_num_1 > I_max_1
       break;
   end
   %%
   signal = H_hat_pre * X;
   %%
   Temp_R = real(signal);
   Temp_I = imag(signal);
   %% Noise Parameter Estimation using Newton's Method
   gamma = func_newton(T_rr,T_lr,T_ri,T_li,gamma_pre,Temp_R,Temp_I);

   iter_num_2 = 0;
   relative_dif_2 = 1;
   H_hat_pre2 = H_hat_pre;
   while(relative_dif_2 > epsilon_2) %% Inner MM iteration
       
       iter_num_2 = iter_num_2 + 1;
       if iter_num_2 > I_max_2
           break;
       end
       signal = H_hat_pre2 * X;
       %%
       Temp_R = real(signal);
       Temp_I = imag(signal);
       
       %% 
       U_rr = T_rr * gamma - Temp_R;                  %% Corresponding to (u - s^hat) in Eq. (13)
       L_rr = T_lr * gamma - Temp_R;                  %% Corresponding to (l - s^hat) in Eq. (13)
       U_ii = T_ri * gamma - Temp_I;                  %% Corresponding to (u - s^hat) in Eq. (13)
       L_ii = T_li * gamma - Temp_I;                  %% Corresponding to (l - s^hat) in Eq. (13)
   
       numerator_r = normpdf(U_rr) - normpdf(L_rr);    %% numerator in Eq. (13)
       numerator_i = normpdf(U_ii) - normpdf(L_ii);    %% numerator in Eq. (13)
   
       denominator_r = fx_cdf(U_rr) - fx_cdf(L_rr);    %% denominator in Eq. (13)
       denominator_i = fx_cdf(U_ii) - fx_cdf(L_ii);    %% denominator in Eq. (13)
       
       Y_hat_r = Temp_R - numerator_r ./ denominator_r;  %% Eq. (17)
       Y_hat_i = Temp_I - numerator_i ./ denominator_i;  %% Eq. (17)
       
       
       Y_hat = Y_hat_r + 1i * Y_hat_i;         %% Eq. (17)
       
       %%
       H_hat = Y_hat * Xt ;                    %% Eq. (20)
       %%
       Dif_2 = H_hat_pre2 - H_hat;
       
       relative_dif_2 = sum(sum( Dif_2.*conj(Dif_2))) / sum(sum( H_hat .* conj(H_hat)));
       H_hat_pre2 = H_hat;    

   end
   Dif_1 = H_hat_pre2- H_hat_pre;
   relative_dif_1 = sum(sum( Dif_1.*conj(Dif_1))) / sum(sum( H_hat_pre.*conj(H_hat_pre)));

   H_hat_pre = H_hat_pre2;
   gamma_pre = gamma;
   
end
H_hat = H_hat_pre / gamma_pre;

h = reshape(H_hat, Nr*Nt, 1);

h = [real(h);imag(h)];



