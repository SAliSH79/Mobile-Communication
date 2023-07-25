%% HW4
% Teacher : Doctor_Mohammadi
%Student-Number : [9723042]
% University: Amirkabir University of Technology

%% Q 4-25

%% Clear recent data
clc;
close all;
clear;
%% Initialization Data
clc;
n = 4; %path loss exponent
sigma = 6 ;
P0 = 0; % 1 mW
d0 = 1; %1 m
Prmin = -118 ; % W
PrHO = -112; %W

D = 1600; %Distance between BS1 and BS2 [m]
d1 = linspace(1,1600,3200); %distance from BS1

%% Forming probablity
clc;
m1 = P0 - 10*n*log10(d1/d0);
m2 = P0 - 10*n*log10((D - d1)/d0);

arg1 = (PrHO - m1)./sigma;
arg2 = (Prmin - m2)./sigma;

Prob1 = 1 - qfunc(arg1);
Prob2 = qfunc(arg2);

Prob_HoF = Prob1.*Prob2 ;

%% Plotting
clc;
figure(1)
plot(d1,Prob_HoF,'. m');hold on;
plot(d1,Prob1,'-- k');hold on;
plot(d1,Prob2,'- b');hold off;
title("Probablity of Handoff at distance d")
axis([1 1600 0 1])
xlabel ('d(m)')
ylabel('Probablity')
grid on
legend('Prob of HandOff','Prob(Pr1<PrHO)','Pr(Pr2>Prmin)','Location','southeast')











