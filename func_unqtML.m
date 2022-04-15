function h = func_unqtML(Y, X)
%% Input
% Y : received signal matrix
% X : transmit signal matrix
%% Output
% h : estimated real-valued channel vector

H = Y * X' / (X * X');

h = H(:);

h = [real(h); imag(h)];