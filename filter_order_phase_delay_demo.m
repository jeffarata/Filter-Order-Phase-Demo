% Jeff Arata
% 12/18/17
%
% This script looks at the effects of filter order on the phase, group
% delay/time delay of a signal.


% Order it, from easiest to hardest:
% 1) group delay (in general, negative derivative of phase shift) vs filter order
% 2) time delay vs filter order
% 3) phase delay (divide phase shift by negative frequency) vs filter order
% 4) phase shift in radians vs filter order (FIR only?)
% 5) filters magnitude response, frequency and phase
% 6) group delay (neg. derivative) vs. frequency (may be constant)
% 7) phase delay (divide by neg frequency) vs. frequency (may be constant)

% Notes to make:
% phase shift = -2*pi*K*f/fs (group delay and frequency can change here)
% time delay = K/fs = (N-1)/(2*fs)
% group delay = phase shift * fs / (2*pi*f)



clear;
clc;

fs = 8000;      % picking a sampling frequency
f = 5;          % pick a frequency for our signal

%% FIR Filters - using fir1

%% Group Delay vs. Filter Order

% In an FIR filter, the phase shift is usually linear in the form -K*w, 
% where K is the group delay and w is frequency. Since group delay is the
% negative frequency derivative of phase shift, we can find K. But it is 
% the average of the number of delayed samples which is:
% 
%   (0 + 1 + 2 + 3 + ... + (N-1) + N )/N = (N(N-1)/2) / N = (N-1)/2
%
% where N is the number of taps in the filter (1 more than its order)
 

max_order = 1000;        % set maximum filter order
order = 1:max_order;    % filter orders
N = order + 1;          % # of taps in filter
K = (N-1)/2;            % group delay

figure(1)
subplot(4,1,1)
plot(order,K)
title('Group Delay of FIR Filter vs. Filter Order')
ylabel('K, # of samples')


%% Time Delay vs. Filter Order

% Time delay is just the group delay, the average number of samples
% delayed, multiplied by the time in between samples, 1/fs

t_delay = K/fs;

figure(1)
subplot(4,1,2)
plot(order,t_delay)
title('Time Delay vs. Filter Order')
ylabel('K/fs, sec.')


%% Phase Shift vs. Filter Order

% Phase shift can be found to be p = -2*pi*K*f/fs
% Here f is a singular frequency

p = -2*pi*K*f/fs;
figure(1)
subplot(4,1,3)
plot(order,p)
title('Phase Shift in Radians vs. Filter Order')
ylabel('p, radians')


%% Phase Delay vs. Filter Order

% Phase delay is the phase shift divided by the negative frequency
%
% Here, since frequency f is constant, we can plot phase delay over filter
% order. Doing so gives us units of seconds for phase delay and in fact
% becomes 2*pi*time_delay

p_delay = p / (-f);

figure(1)
subplot(4,1,4)
plot(order,p_delay)
title('Phase Delay vs. Filter Order')
ylabel('p/(-f), sec.')
xlabel('Filter Order')




%% Simulation of the effects of filtering a signal

% Here we look at the effects of filtering a signal with a successively
% increasing filter order. We'll plot the original signal and the filtered
% signal on the same plot. We'll choose a lowpass filter with a cutoff
% frequency above those in the signal so that we can see the effects on the
% waveform.

% Maybe calculate above features with each filter? Or plot a star or
% diamond at corresponding samples in input and output signal to emphasize
% movement? Use features calculated above to do so.



t = 0:1/fs:1-1/fs;      % plot over a 1 second interval

f = 5;                  % our frequency
x = cos( 2*pi*f*t );    % our signal

%x = x + 0.5*cos(2*pi*9*t);
%x = x + 0.25*cos(2*pi*15*t);

cutoff = fs/4;      % cutoff is half of maximum frequency

point_follow = round(fs/4.5);

for ii = 1:order(end)
    
    b = fir1(ii, cutoff/(fs/2) );       % Filter coefficients
    y = filter(b,1,x);                  % signal output
    y_time = t + t_delay(ii);           % delayed time vector of output
    figure(2)
    plot(t,x,'b', t,y,'r')              % plot input and output of filter
    axis([0 1 -1.5 1.5])
    title('Original Signal (blue) and Filtered Signal (red)')
    hold on
    plot(t(point_follow), x(point_follow), 'o')     % Plot a point to follow
    if mod(ii,2) == 0                               % Plots the actual point
        y_value = y(point_follow+ii/2);             % if even.
        plot(y_time(point_follow), y_value, 'ro')
    elseif (ii>1)                                   % For odd orders, after
        plot(y_time(point_follow), y_value, 'ro')   % first loop, plot the 
    end                                             % previous even y value
    hold off                                        % for less fluctuation.    
    
    % Try plotting a horizontal line between the follow points and indicate
    % the time delay, group delay, or whatever?
    
end

fprintf('We see that phase shift increases linearly (in magnitude) with increasing filter order for an FIR filter.\n')



%{
Notes and findings:

Waves with only odd harmonics can perfectly cancel themselves out with an
appropriately delayed copy of itself. This delay occurs when filtered by an
order that is an odd multiple of fs/f, where f is the fundamental
frequency, or lowest harmonic.

n = 1, 2, 3, 4, 5, ...
order = (2*n-1)*fs/f
group_delay = K = order/2 = ((2*n-1)/2)*fs/f
time_delay = (order/2)*1/fs = (2*n-1)/(2*f)

%}


% IIR Filters

% Bessel Filters



for order = 1:max_order
    
    
    
    
    
end



% Butterworth Filters









% ideas for further investigation

% do the same visualization looping through filters of different orders but
% with an IIR (maybe a chebyshev tye 2).

% try to calculate phase shift and group delay from that?





% Chebyshev Type 1 Filters




% Chebyshev Type 2 Filters






% Elliptic Filters







