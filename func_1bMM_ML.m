function h = func_1bMM_ML(z_bar, X, Nr, t_bar)
%% Input
% z_bar : the real-valued one-bit quantized vector
% X : pilot signal matrix
% Nr : received antenna number
% t_bar : the real-valued threshold vector used
% for one bit quantization
%% Output
% h : estimated real-valued channel vector
%% Reference: [1] F. Liu, X. Shang, and H. Zhu,¡°Effcient majorization-minimization based channel estimation for one-bit massive MIMO systems,¡± 
%% IEEE Transactions on Wireless Communications, vol. 20, no. 6, pp. 3444¨C3457, Jun. 2021.
%% The equation numbers below are corresponding to the equations in [1].
[Nt, K] = size(X); % Nt : transmit antenna number, K : pilot length

t = t_bar(1:K*Nr) + 1i * t_bar(K*Nr+1:2*K*Nr);
T = reshape(t,Nr,K);  %% T = unvec(t)
T_r = real(T);
T_i = imag(T);

z = z_bar(1:K*Nr) + 1i * z_bar(K*Nr+1:2*K*Nr);
Z = reshape(z,Nr,K); %% Z = unvec(z)
Z_r = real(Z);
Z_i = imag(Z);


Xt = X'/(X*X');  %% Eq. (20) (in [1])

Temp = T * Xt;
Temp2 = Temp * X - T;
Denomi=sum(sum(Temp2 .* conj(Temp2))); %% compute the denominator in Eq. (22)

%% Initialization

sigma_hat = 1;            %% initialize the noise parameter estimate
H_hat_pre = zeros(Nr,Nt); %% initialize the channel matrix estimate
% H_hat_pre = randn(Nr,Nt)+1i * randn(Nr,Nt);
epsilon = 1e-5;
relative_dif = 1;

iter_num = 0;
I_max = 50; %% maximum iteration number

while(relative_dif > epsilon)
%     tic
   iter_num = iter_num + 1;
   if iter_num > I_max
       break;
   end
   %%
   signal = H_hat_pre * X;
   %%
   S_r = Z_r .* (real(signal) - T_r) / sigma_hat;  %% Eq.(12a) (matrix version)
   S_i = Z_i .* (imag(signal) - T_i) / sigma_hat;  %% Eq.(12b)
   %%
   U_r = S_r + sqrt(2/pi) ./ erfcx(-S_r / sqrt(2));  %% Eq. (15a)
   U_i = S_i + sqrt(2/pi) ./ erfcx(-S_i / sqrt(2));  %% Eq. (15b)
   %%
   
   Q_r = U_r .* Z_r;
   Q_i = U_i .* Z_i;
   Q = Q_r + 1i * Q_i;   %% Construct Q
   %%
   Qt = Q * Xt;
   A = Qt * X- Q;
   
   gamma_hat = - real(sum(sum(A.*conj(Temp2))))/Denomi; %% Eq. (22)
   
   H_hat = Qt / gamma_hat + Temp;  %% Eq. (19) /gamma
   
   sigma_hat = 1 / gamma_hat;
   
   %%
   Dif = H_hat_pre- H_hat;
   
   relative_dif = sum(sum( Dif.*conj(Dif))) / Nr;

   H_hat_pre = H_hat;
   
end

h = reshape(H_hat, Nr*Nt, 1);
h = [real(h); imag(h)];

