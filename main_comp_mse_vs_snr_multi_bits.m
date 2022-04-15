clear all;
close all;
clc;

Nr = 64;                    %% receive antenna number
Nt = 32;                    %% transmit antenna number
K = 64;                     %% pilot length
SNR_dB = -10:5:15;          %% SNR
SNR = 10.^(SNR_dB / 10);

Pnt = 1 ./ SNR;                 %% noise power (the transmit power is normalized as 1)
MonteC = 10;                    %% Monte Carlo trail number

%% Rayleigh Fading Channel Model

H = randn(Nr,Nt) + 1i* randn(Nr,Nt);
H = H/sqrt(2);
h_true = H(:);
h_true = [real(h_true); imag(h_true)];
H_fro = norm(H,'fro')^2;


%% MSE vector for different quantizaton bit number

mse_1bit_mc = zeros(length(SNR), MonteC);
mse_2bit_mc = zeros(length(SNR), MonteC);
mse_3bit_mc = zeros(length(SNR), MonteC);
mse_4bit_mc = zeros(length(SNR), MonteC);
mse_unqt_mc = zeros(length(SNR), MonteC); 

%% Zadoff-Zhu sequence / Orthogonal Pilots
x = zadoff_chu(K);
S0 = zeros(K,K);
ind = 1:K;
for jj = 1:K
    ind = mod(ind,K) + 1;
    S0(:,jj) = x(ind);
end

for i = 1 : length(SNR)
    Pn = Pnt(i);
    
    %% Thresholds used for 1-bit quantization
    h_max = sqrt(1 + Pn)/sqrt(2);
    ht = -h_max : 2*h_max/7 : h_max;
    ht = ht.';
    
    for j = 1: MonteC
%         [i, j]
        %% Thresholds used for 1-bit quantization
        t_bar = zeros(2*K*Nr,1);
        kro = ones(K,1);
        tem = ht(randi(8, Nr, 1));
        t_bar(1:K*Nr) = kron(kro, tem);
        tem = ht(randi(8,Nr,1));
        t_bar(K*Nr+1:2*K*Nr) = kron(kro, tem);
        
        %% Received Signal Model
        %   qpsk_signal = randi(4,Nt,K);
        %   X = pskmod(qpsk_signal - 1, 4, pi/4);        %% Random QPSK pilots
        X = S0(randperm(Nt),:);                          %% Random Orthogonal Pilots
        X = X / sqrt(Nt);                                %% normalize the transmit power to 1
        N = random('norm',0,sqrt(Pn/2),Nr,K) + 1i * random('norm',0,sqrt(Pn/2),Nr,K);  %% Noise matrix
        Y_unqt = H * X  + N;
        y_unqt = Y_unqt(:);
        y_bar = [real(y_unqt); imag(y_unqt)];
        
        %% 1-bit
        z_bar = sign(y_bar - t_bar);
        y_lowb = z_bar(1 : K*Nr) + 1i * z_bar(K*Nr+1 : 2*K*Nr);
        Y_lowb = reshape(y_lowb, Nr, K); 
        h_1bit = func_1bMM_ML(z_bar, X, Nr, t_bar);
        mse_1bit_mc(i, j) = sum((h_1bit - h_true).^2);
        
        
        %% 2-bit
        [y_lowb, thres_r, thres_l] = func_lowbSamp(y_bar, 2);
        h_2bit = func_CQML(thres_l, thres_r, y_lowb, X, Nr);
        mse_2bit_mc(i, j) = sum((h_2bit - h_true).^2);
        
        
        %% 3-bit
        [y_lowb, thres_r, thres_l] = func_lowbSamp(y_bar, 3);
        h_3bit = func_CQML(thres_l, thres_r, y_lowb, X, Nr);
        mse_3bit_mc(i, j) = sum((h_3bit - h_true).^2);
        
        %% 4-bit
        [y_lowb, thres_r, thres_l] = func_lowbSamp(y_bar, 4);
        h_4bit = func_CQML(thres_l, thres_r, y_lowb, X, Nr);
        mse_4bit_mc(i, j) = sum((h_4bit - h_true).^2);
        
        %% Unquantized case
        h_unqt = func_unqtML(Y_unqt, X);
        mse_unqt_mc(i, j) = sum((h_unqt - h_true).^2);
        
    end
    
end

%% Compute the NMSE
mse_1bit = mean(mse_1bit_mc, 2)/H_fro;
mse_2bit = mean(mse_2bit_mc, 2)/H_fro;
mse_3bit = mean(mse_3bit_mc, 2)/H_fro;
mse_4bit = mean(mse_4bit_mc, 2)/H_fro;
mse_unqt = mean(mse_unqt_mc, 2)/H_fro;

%% Plot figure
fx=@(x) 10*log10(x);
A=[0.13 0.55 0.13;0.60 0.20 0.98]; 
figure(1);
hold on;
plot(SNR_dB,fx(mse_1bit),'-*','COLOR',A(1,:),'LineWidth',1.5,'MarkerSize',14);
plot(SNR_dB,fx(mse_2bit),'-ks','LineWidth',1.5,'MarkerSize',14);
plot(SNR_dB,fx(mse_3bit),'-kd','LineWidth',1.5,'MarkerSize',14);
plot(SNR_dB,fx(mse_4bit),'-ko','LineWidth',1.5,'MarkerSize',14);
plot(SNR_dB,fx(mse_unqt),'-r>','LineWidth',1.5,'MarkerSize',14);
grid on; 
set(gca, 'FontSize', 14);
lg = legend('1bMM-ML','CQML:2-bit','CQML:3-bit','CQML:4-bit','UnqtML');
lg.FontSize = 14;
xlabel('SNR (dB)','FontSize',14);
ylabel('MSE (dB)','FontSize',14);

