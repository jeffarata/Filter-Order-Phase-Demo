% Jeff Arata
% 1/5/18

function [ D, W ] = group_delay( b, a, N, fs )
% This function calculates the group delay of a filter based on its
% coefficients using the method laid out in Julius Smith's book Intro to
% Digital Filters under the section on Numerical Computation of Group
% Delay, returning group delay in samples.

% Inputs:
%
% b -       feedforward coefficients
% a -       feedback coefficients
% n -       number of data points/length of return values, D and W
% fs -      sampling frequency, Hz
%
% Outputs:
% 
% D -       Group delay in samples
% W -       angular frequency vector, returns as true frequency (Hz) if fs
%           is specified

if nargin < 2
    a = 1;              % Default to an FIR filter
    N = 512;            % Default number of data points for output
    fs = 2;             % Allows calculation of angular frequency below
elseif nargin < 3
    N = 512;            % Default number of data points for output
    fs = 2;             % Allows calculation of angular frequency below
elseif nargin < 4
    fs = 2;             % Allows calculation of angular frequency below
else
end

fft_size = 2*(N-1);

B_len = length(b);
A_len = length(a);

Br = fft(b.*[0:B_len-1], fft_size);         % fft of ramped coefficients b
B = fft(b, fft_size);
Ar = fft(a.*[0:A_len-1], fft_size);         % fft of ramped coefficients a
A = fft(a, fft_size);

D =  real( Br./B ) - real( Ar./A );         % Group delay in samples
D = D(1:fft_size/2+1);                      % Group delay from 0 to pi

W = linspace(0, 2*pi*((N-1)/N)*(fs/2), N);  % angular frequency vector
if nargin == 4
    W = linspace(0, ((N-1)/N)*(fs/2), N);   % True frequency (Hz) if fs specified
end

end
