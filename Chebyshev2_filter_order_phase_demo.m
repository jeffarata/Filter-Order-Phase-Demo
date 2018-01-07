% Jeff Arata
% 1/2/18
%
% This script looks at the effects of filter order on the phase, group
% delay/time delay, and phase delay of a signal, using a Chebyshev Type 2 
% lowpass filter.


% Notes on IIR filters:
%
% phase response = angle(H); H is complex frequency response
% group delay in samples = -d(phase)/d(w); w is angular frequency
% group delay in words: negative derivative of phase shift
% phase delay = -phase ./ w = -phase ./ (2*pi*f)
% phase delay in words: negative phase shift divided by angular frequency

close all;
clear;
clc;


%% Simulation of the effects of the Chebyshev Type 2 Lowpass Filter. 

% Here we look at the effects of filtering a signal with a successively
% increasing filter order. We'll plot the original signal and the filtered
% signal on the same plot. We'll choose a lowpass filter with a cutoff
% frequency above those in the signal so that we can see the effects on the
% waveform.

fs = 10000;             % picking a sampling frequency
f = 5;                  % pick a frequency for our signal

t = 0:1/fs:1-1/fs;      % plot over a 1 second interval
x = cos( 2*pi*f*t );    % our signal

IIR_max_order = 25;
cutoff = (fs/4);                 % choose filter cutoff frequency <= fs/2
N = 2048;                        % number of data points wanted
cutoff_idx = round((cutoff/(fs/2))*N + 1);

idx_check = round(linspace(1,cutoff_idx,5));
delay_check_matrix = zeros(IIR_max_order,length(idx_check));
phase_check_matrix = zeros(IIR_max_order,length(idx_check));


%% Chebyshev Type 2 Filters

R = 20;         % number of decibels down the stopband ripple is

for ii = 1:IIR_max_order
    
    [b,a] = cheby2(ii, R, cutoff/(fs/2) );  % Filter coefficients
    y = filter(b,a,x);                      % signal output
    % 1) Signal and filtered output
    figure(1)
    plot(t,x,'b', t,y,'r')              % plot input and output of filter
    axis([0 t(end)+1/fs -1.5 1.5])
    title('Original Signal (blue) and Chebyshev2 Filtered Signal (red)')    
    
    % 2) Frequency Response of filter
    figure(2)
    [H,W] = freq_response(b, a, N, fs); % W is true frequency, Hz
    subplot(2,1,1)
    title(['Frequency Response Chebyshev2 Order: ', num2str(ii)])
    subplot(2,1,2)
    axis([0 fs/2 -3 0])       
    
    % Returns group delay in samples and W in Hz
    [D, D_W] = group_delay(b, a, N, fs);
    D = D/fs;  % Group delay in seconds
    % Phase response of the filter
    phase = angle(H);
    phase = unwrap(phase);                  
            
    % 3.1) Group Delay in Seconds
    figure(3)
    subplot(3,1,1)
    plot(D_W, D)
    ylabel('Group Delay, sec')
    title(['Group Delay vs. Frequency - Chebyshev2 Order: ', num2str(ii)])
    axis([0 fs/2 0 0.03])
    
    % 3.2) Group Delay in Samples
    sample_delay = D / (1/fs);  % Number of samples delayed per freq.
    figure(3)
    subplot(3,1,2)
    plot(D_W, sample_delay)
    ylabel('# of Samples Delay')
    title(['Samples Delayed vs. Frequency - Chebyshev2 Order: ', num2str(ii)])
    axis([0 fs/2 0 300])
    
    % 3.3) Phase Delay in Seconds
    phase_delay = -phase ./ (W*2*pi);   % time delay of phase
    figure(3)                           % multiply by 2*pi for angular frequency
    subplot(3,1,3)
    plot(W, phase_delay)
    xlabel('Frequency, Hz')
    ylabel('Phase Delay, sec')
    title(['Phase Delay vs. Frequency - Chebyshev2 Order: ', num2str(ii)])
    axis([0 fs/2 0 0.0006])
        
    % For plotting the delay of a few frequencies vs. filter order
    delay_check_matrix(ii,:) = D(idx_check);
    % For plotting the phase shift of a few frequencies vs. filter order
    phase_check_matrix(ii,:) = phase(idx_check);
    
end

% 4.1) Plot delay of certain frequencies vs filter order
figure(4)
subplot(2,1,1)
plot(delay_check_matrix)
title(['Delay in Secs vs. Chebyshev2 Order of ', num2str(length(idx_check)), ' Different Frequencies'])
xlabel('Filter Order')
ylabel('Delay in Seconds')
freq_legend = W(idx_check);
legend(num2str(freq_legend', '%5.0f\n'), 'Location', 'NorthWest')
legend BOXOFF

% 4.2) Plot the change in delay of certain frequencies vs filter order
diff_delay_check_matrix = diff(delay_check_matrix);
figure(4)
subplot(2,1,2)
plot(diff_delay_check_matrix)
title(['Diff of Delay in Secs at ', num2str(length(idx_check)), ' Frequencies vs. Chebyshev2 Order'])
xlabel('Filter Order')
ylabel('Difference of Delay in Seconds')

% Note on Group Delay:
%
% From these plots, we see that the difference in group delay for
% frequencies in the passband approaches 0 with increasing filter order 
% (past order 5 or so). In fact, the difference approaches 0 from below, so
% the group delay starts out increasing and then decreases towards 0. As 
% for the cutoff frequency, the difference in delay is increasing, as we 
% would expect it to. For low filter order values (<5), the difference in
% delay time is erratic.

% 5.1) Plot the phase shift vs filter order of certain frequencies
figure(5)
subplot(2,1,1)
plot(phase_check_matrix/pi)
title(['Phase Shift (pi radians) vs. Chebyshev2 Order of ', num2str(length(idx_check)), ' Different Frequencies'])
xlabel('Filter Order')
ylabel('Phase Shift (pi radians)')
freq_legend = W(idx_check);
legend(num2str(freq_legend', '%5.0f\n'), 'Location', 'NorthWest')
legend BOXOFF

% 5.2) Plot the change in phase shift of certain frequencies vs filter order
diff_phase_check_matrix = diff(phase_check_matrix);
figure(5)
subplot(2,1,2)
plot(diff_phase_check_matrix/pi)
title(['Diff of Phase Shift (pi radians) at ', num2str(length(idx_check)), ' Frequencies vs. Chebyshev2 Order'])
xlabel('Filter Order')
ylabel('Difference of Phase Shift (pi radians)')

% Note on Phase Shift:
% 
% For frequencies in the passband, we see the difference in phase shift
% with increasing filter order approach 0. For frequency bin 0, it is
% constant at 0. For others, the phase shift seems to start out by delaying
% more and more with filter order and then turns around back towards 0. For
% frequencies in the passband, the phase shift approaches a constant,
% especially after filter order 10 or so.

% Note on the Cutoff Frequency:
%
% We see the difference in phase shift approach 0 for the cutoff frequency.
% It still remains a delay such that the phase shift is always increasingly
% negative. So the delay with increasing filter order at the cutoff ranges
% from 3pi/10 to approaching 0 radians.
