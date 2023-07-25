%% HW4

% University: Amirkabir University of Technology

%% Q 4-31
%% Clear recent data
clc;
close all;
clear;
%% Initialization
dmax = 3000; %as Example 4.9
n = 4.4; %as Example 4.9
sigma = 6.17; %As Example 4.9 in db
d0 = 100;%as Example 4.9
PL_d0 = -20;%as Example 4.9
n_sample = 750;%as Example 4.9
%% MSE
if dmax >= d0
        % caculate the path loss
        d = linspace(d0, dmax, n_sample);
        X = normrnd(0, sigma, [1, n_sample]);
        PL = PL_d0 + 10*n*log10(d/d0) + X;
        
        % plotting the results
        scatter(PL,20*log(d),"r");
        title("Path Loss / distance")
        xlabel ('d(m)')
        ylabel('Path Loss Mag in db')
        grid on
        legend('Path Loss(db)','Location','southeast')
end






