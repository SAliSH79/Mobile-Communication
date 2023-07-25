%% HW4
% University: Amirkabir University of Technology

%% Q 5-20
%% Clear recent data
clc;
close all;
clear;
%% Initialization
fm = 200;
p = linspace(0,1,1000);
rayChan = comm.RayleighChannel('SampleRate',10000,'MaximumDopplerShift',fm);
sig = j*ones(2000,1); % Signal
out = rayChan(sig); % Pass signal through channel.
%rayChan % Display all properties of the channel object.

%% Plot power of faded signal, versus sample number.
figure(1)
plot(20*log10(abs(out)))
title("Simulated Path Loss Faded")
xlabel ('d(m)')
ylabel('Path Loss Mag in db')
grid on
legend('Faded Path Loss --> fd = 200Hz','Location','southeast')

%% Calculation
Nr = sqrt(2*pi)*fm.*(p).*exp(-p.^2); %level Crossing rate
taw = (exp(p.^2) - 1)./(p*fm*sqrt(2*pi)) ; %average fade durations

%% Plotting
clc;
figure(2)


subplot(211)
plot(p,Nr,'. m');
title("level Crossing rates")
xlabel ('p')
ylabel('Nr(s)')
grid on
legend('Crossing rate','Location','southeast')

subplot(212)
plot(p,taw,'-- k');
title("Average Fade")
xlabel ('p')
ylabel('taw(s)')
grid on
legend('Average Fading','Location','southeast')

