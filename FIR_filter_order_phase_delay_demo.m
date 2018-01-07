% Jeff Arata
% 12/18/17
%
% This script looks at the effects of filter order on the phase, group
% delay/time delay, and phase delay of a signal, using an FIR filter.


% Notes on FIR Filters
% phase shift = -2*pi*K*f/fs (group delay and frequency can change here)
% group delay in samples = K = (N-1)/2 = phase shift * fs / (2*pi*f)
% group delay in time = K/fs = (N-1)/(2*fs)
% phase delay = phase shift / (2*pi*f); Phase divided by angular frequency

close all;
clear;
clc;


%% FIR Filters - using fir1

fs = 10000;      % picking a sampling frequency
f = 10;          % pick a frequency for our signal


%% 1.1) Group Delay vs. Filter Order

% In an FIR filter, the phase shift is usually linear in the form -K*w, 
% where K is the group delay and w is frequency. Since group delay is the
% negative frequency derivative of phase shift, we can find K. But it is 
% the average of the number of delayed samples which is:
% 
%   (0 + 1 + 2 + 3 + ... + (N-1) + N )/N = (N(N-1)/2) / N = (N-1)/2
%
% where N is the number of taps in the filter (1 more than its order)
 

FIR_max_order = 100;        % set maximum filter order
order = 1:FIR_max_order;    % filter orders
N = order + 1;          % # of taps in filter
K = (N-1)/2;            % group delay

figure(1)
subplot(4,1,1)
plot(order,K)
title('Group Delay (Samples) of FIR Filter vs. Filter Order')
ylabel('K, # of samples')


%% 1.2) Time Delay vs. Filter Order

% Time delay is just the group delay, the average number of samples
% delayed, multiplied by the time in between samples, 1/fs

t_delay = K/fs;

figure(1)
subplot(4,1,2)
plot(order,t_delay)
title('Group Delay (Time) vs. Filter Order')
ylabel('K/fs, sec.')


%% 1.3) Phase Shift vs. Filter Order

% Phase shift can be found to be p = -2*pi*K*f/fs
% Divide out pi for units of pi radians
% Here f is a singular frequency in Hz, not angular frequency
   
p = -2*pi*K*f/fs;       % 2*pi*f is angular frequency

figure(1)
subplot(4,1,3)
plot(order,p/pi)
title(['Phase Shift in Pi Radians vs. Filter Order for Frequency = ', num2str(f), ' Hz'])
ylabel('p, pi radians')


%% 1.4) Phase Delay vs. Filter Order

% Phase delay is the phase shift divided by the negative angular frequency
%
% Here, since frequency f is constant, we can plot phase delay over filter
% order. Doing so gives us units of seconds for phase delay and in fact
% becomes 2*pi*time_delay

p_delay = p / (-f*2*pi);    % 2*pi*f is angular frequency

figure(1)
subplot(4,1,4)
plot(order,p_delay)
title(['Phase Delay vs. Filter Order for Frequency = ', num2str(f), ' Hz'])
ylabel('p/(2*pi*(-f)), sec.')
xlabel('Filter Order')


%% Simulation of the effects of filtering a signal

% Here we look at the effects of filtering a signal with a successively
% increasing filter order. We'll plot the original signal and the filtered
% signal on the same plot. We'll choose a lowpass filter with a cutoff
% frequency above those in the signal so that we can see the effects on the
% waveform.

t = 0:1/fs:1-1/fs;      % plot over a 1 second interval
x = cos( 2*pi*f*t );    % our signal

N_data = 2048;      % Number of data points for output 
cutoff = fs/4;      % cutoff is half of maximum frequency fs/2
cutoff_idx = round((cutoff/(fs/2))*N_data+1);

% initialize for seeing how delay and phase changes with filter order for
% certain frequencies
idx_check = round(linspace(1,cutoff_idx,5));
delay_check_matrix = zeros(FIR_max_order,length(idx_check));
phase_check_matrix = zeros(FIR_max_order,length(idx_check));

% Index to plot a point on input and output signal to help visualize delay
point_follow = round(fs/4.5);   

for ii = 1:order(end)
    
    b = fir1(ii, cutoff/(fs/2) );       % Filter coefficients
    y = filter(b,1,x);                  % signal output
    y_time = t + t_delay(ii);           % delayed time vector of output
    % 2) Plot of input and output signals
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
    
    % 3) Frequency Response
    figure(3)
    [H, W] = freq_response(b, 1, N_data, fs);     % W is true frequency, Hz
    subplot(2,1,1)
    title(['Frequency Response FIR Order: ', num2str(ii)])
    
    % Returns group delay in samples and W in Hz
    [D, D_W] = group_delay(b, 1, N_data, fs);
    D = D/fs;   % Group delay in seconds
    % For plotting the delay of a few frequencies vs. filter order
    delay_check_matrix(ii,:) = D(idx_check);

    % Phase of filter
    phase = angle(H);
    phase = unwrap(phase);
    % For plotting the phase shift of a few frequencies vs. filter order
    phase_check_matrix(ii,:) = phase(idx_check);
    
end

fprintf('We see that phase shift increases linearly (in magnitude) with increasing filter order for an FIR filter.\n')

% 4.1) Plot delay of certain frequencies vs filter order
figure(4)
subplot(2,1,1)
plot(delay_check_matrix);
title(['Delay in Seconds vs. FIR1 Order of ', num2str(length(idx_check)), ' Different Frequencies'])
xlabel('Filter Order')
ylabel('Delay in Seconds')
freq_legend = W(idx_check);
legend(num2str(freq_legend', '%5.0f\n'), 'Location', 'NorthWest')
legend BOXOFF

% 4.2) Plot the change in delay of certain frequencies vs filter order
diff_delay_check_matrix = diff(delay_check_matrix);
figure(4)
subplot(2,1,2)
plot(diff_delay_check_matrix);
title(['Diff of Delay in Seconds at ', num2str(length(idx_check)), ' Frequencies vs. FIR1 Order'])
xlabel('Filter Order')
ylabel('Difference of Delay in Seconds')

% Notes on Group Delay:
%
% From these graphs, we see that the difference in delay time of an
% increasing filter order is the same constant value across frequencies in
% the passband. This we would expect due to the linear phase of an FIR.
% There is quite a bit of noise with increasing filter order, but this can
% be attributed to numerical error be observing the y-axis range.


% 5.1) Plot the phase shift vs filter order of certain frequencies
figure(5)
subplot(2,1,1)
plot(phase_check_matrix/pi)     % Divide out pi for units of pi radians
title(['Phase Shift (pi radians) vs. FIR1 Order of ', num2str(length(idx_check)), ' Different Frequencies'])
xlabel('Filter Order')
ylabel('Phase Shift (pi radians)')
freq_legend = W(idx_check);
legend(num2str(freq_legend', '%5.0f\n'), 'Location', 'NorthWest')
legend BOXOFF

% 5.2) Plot the change in phase shift of certain frequencies vs filter order
diff_phase_check_matrix = diff(phase_check_matrix);
figure(5)
subplot(2,1,2)
plot(diff_phase_check_matrix/pi)    % Divide out pi for units of pi radians
title(['Diff of Phase Shift (pi radians) at ', num2str(length(idx_check)), ' Frequencies vs. FIR1 Order'])
xlabel('Filter Order')
ylabel('Difference of Phase Shift (pi radians)')

% Note on Phase Shift:
%
% Observing the difference in phase shift vs increasing filter order, we
% see that this is constant for frequencies in the passband of an FIR
% filter. This makes sense with an FIR's linear phase. In addition, the
% value of the phase shift for frequencies seems to be spread linearly,
% which also makes sense.


% Note on the Cutoff Frequency:
% 
% At the cutoff frequency, we see that increasing the filter order by 1
% always increases the phase shift by pi/4 radians!


%{
Other Notes and Findings:

Waves with only odd harmonics can perfectly cancel themselves out with an
appropriately delayed copy of itself. This delay occurs when filtered by an
order that is an odd multiple of fs/f, where f is the fundamental
frequency, or lowest harmonic.

n = 1, 2, 3, 4, 5, ...
order = (2*n-1)*fs/f
group_delay = K = order/2 = ((2*n-1)/2)*fs/f
time_delay = (order/2)*1/fs = (2*n-1)/(2*f)
%}
