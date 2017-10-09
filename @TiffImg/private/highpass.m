function x_f = highpass(x,T,dx)
% HIGHPASS Apply high pass filter to signal x by removing all period
% components less than T
%
% x_f = highpass(x,T)
%
% Input 
%   x : signal
%   T : period cut-off for high pass filter
%   dx : sample period (1/dx is the sample frequency)
%
% Output
%  xh : filtered signal

% James Kapaldo


% Cut-off frequency index.
L = size(x,1);
freq_cut_off_idx = round(dx*L/T);

% Frequency domain
xh = fft(x);

% Remove frequecies lower than the cut-off
xh([1:freq_cut_off_idx, end-freq_cut_off_idx+1:end],:) = 0;

% Back to spatial space.
x_f = real(ifft(xh));


