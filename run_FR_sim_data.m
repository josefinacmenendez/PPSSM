% run_FR_sim_data.m
% Written by Josefina Correa
% This script loads simulated data from 'sim_data.mat', and then computes
% spike rates via state space smoothing
% sim_data.mat contains:
% Y: a binary matrix of a time-series of 1000 simulated spike trains, of
% 100 second in duration
% lambda: the true state underlying the observations
close all;
clear all;
load('sim_data.mat');
I = sum(Y); % dims 1 x Time
J = size(Y,1);
Fs = 1000; 
initial_guess = 1;
sigesqguess = 1e-3;
sigesqkguess = 1e-3;
xkguess = log(mean(I));
% estimate starting guess of first time-point
FR_initial_guess = get_FR_function(I, J,initial_guess, sigesqguess, sigesqkguess, xkguess);
% use initial guess to estimate all other time-points
FR_struct = get_FR_function(I,J, false, FR_initial_guess.sige_sq, FR_initial_guess.sig_sq_0_K, FR_initial_guess.x_0_K);
% Compute confidence intervals
N_draws = 1000; %1000 Monte Carlo samples
alpha = 0.05; %confidence level
fr_function_cis = compute_fr_function_cis(FR_struct, N_draws, alpha);

T = [2:numel(I)]./Fs;
x_k_K = exp(FR_struct.x_k_given_K) .*Fs;
FR_cis= exp(fr_function_cis) .* Fs;
state_per_second = lambda(2:end) .* Fs;
emp_rate_per_sec= I(2:end)./J .* Fs;
f = figure('renderer','painters');

plot(T, x_k_K, '--k','linewidth', 1.5); hold on;
plot(T, state_per_second, '--r','linewidth', 1.5); hold on;
plot(T, emp_rate_per_sec,'--', 'color', [0 0 1 0.15],'linewidth', 1.5); hold on;
plot(T, FR_cis(1,:), '--', 'color', 'b','linewidth', 1.5); hold on;
plot(T, FR_cis(3,:),'--', 'color', 'b','linewidth', 1.5);
legend('Estimated State', 'True State', 'Empirical rate' ,'95% CI' , 'location', 'southoutside', 'orientation', 'horizontal');
xlabel('Time (seconds)');
ylabel('Spikes per second');
axs = gca;
axs.FontSize = 20;

