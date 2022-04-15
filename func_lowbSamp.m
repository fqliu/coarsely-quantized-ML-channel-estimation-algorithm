function [y,h,l] = func_lowbSamp(s, L)
%% L bits Lloyd-Max quantizer 
%% Input
% s : input signal waitting for sampling (real-valued)
% L : number of sampling bits
%% Output
% y : L-bit quantized sigal
% h : upper bound of samples
% l : lowwer bound of samples
%% Reference: J. Max, ¡°Quantizing for minimum distortion,¡± IRE Trans. Inf. Theory, vol. 6, no. 1, pp. 7¨C12, Mar. 1960.
ampr = sqrt(mean(s.^2));
s = s/ampr;
switch L
    case 1, % 2=2^1 bins
        thr = [-Inf, 0, Inf];
        lev = [-0.7979, 0.7979];
    case 2, % 4=2^2 bins
        thr = [-1e10   -0.9816    0.0000    0.9816       1e10];
        lev = [-1.5104   -0.4528    0.4528    1.5104];
    case 3, % 8=2^3 bins
        thr = [-1e10   -1.7479   -1.0500   -0.5005   -0.0000 ...
            0.5005    1.0500    1.7479       1e10];
        lev = [-2.1519   -1.3439   -0.7560   -0.2451    0.2451 ...
            0.7560    1.3439    2.1519];
    case 4, % 16=2^4 bins
    thr = [-1e10   -2.401   -1.844   -1.437   -1.099  ...
            -0.7996   -0.5224   -0.2582    0.0000    0.2582  ...
             0.5224    0.7996    1.099    1.437    1.844  ...
             2.401       1e10];
         lev = [-2.733   -2.0690   -1.6180   -1.256   -0.9424  ...
             -0.6568   -0.3881   -0.1284    0.1284    0.3881    0.6568 ...
             0.9424    1.256   1.6180    2.0690    2.733];
    otherwise
        error('This demo file works only for 1 <= L <= 4');
end

thr_begin = thr(1:end-1);
thr_end   = thr(2:end);   

Qc = @(x) lev([x >= thr_begin] & [x < thr_end]); 
Qh = @(q) thr_end([q == lev(1:end)]);
Ql = @(q) thr_begin([q == lev(1:end)]);

Q     = @(x) arrayfun(Qc,x);
thr_h = @(q) arrayfun(Qh,q);
thr_l = @(q) arrayfun(Ql,q);

yr0 = Q(s);      yr = yr0*ampr;
hr0 = thr_h(yr0);  hr = hr0*ampr;
lr0 = thr_l(yr0);  lr = lr0*ampr;
y=yr;
h=hr;
l=lr;

