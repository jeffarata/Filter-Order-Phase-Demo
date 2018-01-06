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
% a     - feedback coefficients
% N     - # of data points for output vectors - half fft size
% fs    - the sampling frequency
%
% Outputs:
%
% H     - the frequency response, complex valued
% W     - the angular frequency vector; if fs is specified, W is true
%         frequency

if nargin < 2
    a = 1;          % defaults a single feedback coefficient
    N = 512;        % default # of data point for output
    fs = 2;         % allows unaffected normalized frequency vector
    xlabel_string = 'Normalized Frequency';
elseif nargin < 3
    N = 512;        % defualt # of data point for output
    fs = 2;         % allows unaffected normalized frequency vector
    xlabel_string = 'Normalized Frequency';
elseif nargin < 4
    fs = 2;         % allows unaffected normalized frequency vector  
    xlabel_string = 'Normalized Frequency';
else
    xlabel_string = 'True Frequency, Hz';  
end
    
fft_size = 2*(N-1);
   
H_full = fft(b, fft_size);            % freq response for FIR case
if length(a) > 1
    H_full = H_full ./ fft(a, fft_size); % for IIR case
end
H = H_full(1:fft_size/2+1);            % frequency response at freq. [0, pi]

W = linspace(0, 2*pi*((N-1)/N)*(fs/2), N);  % angular frequency vector
if nargin == 4
    W = linspace(0, ((N-1)/N)*(fs/2), N);   % True frequency (Hz) if fs specified
end

% Plotting
magH = abs(H);
magH = 20*log10(magH/max(magH));       % Magnitude in dB
phase = angle(H);
phase = unwrap(phase)/pi;              % unwrapped phase in pi*radians
      
subplot(2,1,1)
plot(W, magH)
ylabel('Magnitude Response, dB')
xlabel(xlabel_string)
axis([W(1) W(end) -120 10])
subplot(2,1,2)
plot(W, phase)
ylabel('Phase in Pi Radians')
xlabel(xlabel_string)
axis([W(1) W(end) (min(phase)-0.15*abs(min(phase))) 0])

end
