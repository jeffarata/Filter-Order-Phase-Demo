% Jeff Arata
% 12/25/17
% This function finds the frequency response of a filter, FIR or IIR, plots
% the magnitude and phase of the response and returns the frequency
% response as a vector along with the frequency vector.

function [ H, W ] = freq_response( b, a, N, fs )
% Return and plot the frequency response of a filter.
%
% Inputs:
%
% b     - feedforward coefficients
% a     - foodback coefficients
% N     - number of desired data points for the frequency reponse
%
% Outputs:
%
% H     - the frequency response, complex valued
% W     - the frequency vector

if nargin < 2
    a = 1;          % defaults a single feedback coefficient
    N = 512;        % default size of FFT
    fs = 2;         % allows unaffected normalized frequency vector
    xlabel_string = 'Normalized Frequency';
    
elseif nargin < 3
    N = 512;        % defualt size of FFT
    fs = 2;         % allows unaffected normalized frequency vector
    xlabel_string = 'Normalized Frequency';
    
elseif nargin < 4
    fs = 2;         % allows unaffected normalized frequency vector  
    xlabel_string = 'Normalized Frequency';
    
else
    xlabel_string = 'True Frequency, Hz';
    
end
    
b = [b, zeros(1, N-length(b))];    
H_full = fft(b);

if length(a) > 1
    a = [a, zeros(1, N-length(a))];
    H_full = H_full ./ fft(a);
end

H = H_full(1:N/2+1);            % frequency response at freq. [0, pi]

magH = abs(H);
magH = 20*log10(magH/max(magH));    % Magnitude in dB
phase = angle(H);
phase = unwrap(phase);              % unwrapped phase in radians

W = (0:N/2).*(fs/2)./(N/2);         % frequency vector

subplot(2,1,1)
plot(W, magH)
ylabel('Magnitude Response, dB')
xlabel(xlabel_string)
axis([W(1) W(end) -150 10])
subplot(2,1,2)
plot(W, phase)
ylabel('Phase in Radians')
xlabel(xlabel_string)

end

