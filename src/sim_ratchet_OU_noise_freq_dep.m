%SIM_RATCHET_OU_NOISE_FREQ_DEP simulates information ratchet with additonal
%nonequilibrium Ornstein-Uhlenbeck noise and plots output power as a
%function of correlation frequency
%
% OUTPUTS:
%  - creates eps figures of output power as function of noise correlation
%    frequency
%  - creates csv-file of correlation frequency vs. output power
%
% author:  JEhrich
% version: 1.4 (2022-07-28)
% changes: updated experimental parameters
clear
close all
clc
% set font size, line width, and marker size
fS = 18;
lW = 2.5;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% system parameters
Dne = 3.04;
fne_vec = logspace(-3,6,1E2+1);
% sampling time
ts = 1/39;
% scaled effective mass
dg = 0.37;

%% simulation parameters
% number of samples
K = 1E4;
% number of runs
N = 1E4;

%% main loop
% data structures for mean and variance of velocity
F_dot = nan(length(fne_vec),1);
var_F_dot = nan(length(fne_vec),1);
tic
parfor ii = 1:length(fne_vec)
    ii
    fne = fne_vec(ii);
    % one simulation to equilibrate
    [x_traj,l_traj,zeta_traj] = sim_OU_ratchet(dg,Dne,fne,ts,K,2);
    % run N simulations, save final positions
    x_vec = nan(N,1);
    l_vec = nan(N,1);
    for jj = 1:N
        % re-center around l = 0
        x = x_traj(end) - l_traj(end);
        l = 0;
        [x_traj,l_traj,zeta_traj] = sim_OU_ratchet(dg,Dne,fne,ts,K,2,x,l,zeta_traj(end));
        % write out positions
        x_vec(jj) = x_traj(end);
        l_vec(jj) = l_traj(end);
    end
    % compute output power, mean and variance
    F_dot(ii) = mean(l_vec*dg/K/ts);
    var_F_dot(ii) = var(l_vec*dg/K/ts);
    
end
toc

%% compute white-noise limit power
F_dot_white = dg*sqrt(2)*sqrt(1 + Dne)*exp(-dg^2/(2 + 2*Dne))/(sqrt(pi)*(1 + erf(dg*sqrt(2)/(2*sqrt(1 + Dne)))));

%% power without active noise
F_dot_0 = dg*sqrt(2)*exp(-dg^2/2)/(sqrt(pi)*(1 + erf(sqrt(2)*dg/2)));

%% plot output power and simulation data
figure();
errorbar(fne_vec,F_dot,sqrt(var_F_dot/N),'bs','Linewidth',lW,'MarkerSize',mS);
hold on;
plot(fne_vec,F_dot_white*ones(size(fne_vec)),'k:','Linewidth',lW,'MarkerSize',mS);
plot(fne_vec,F_dot_0*ones(size(fne_vec)),'k--','Linewidth',lW,'MarkerSize',mS);
xlabel('O-U noise correlation frequency $f_\mathrm{neq}$','interpreter','latex');
ylabel('rate of free energy gain','interpreter','latex');
set(gca,'XScale','log','FontSize',fS);
title(['sampling time $t_\mathrm{s}=' num2str(ts) '$, noise strength $D_\mathrm{ne}=' num2str(Dne) '$' ],'FontWeight','Normal','interpreter','latex');
legend({'simulation','white-noise limit','no active noise'},'Location','NorthWest');
% save figure
saveas(gcf, '../doc/ratchet_OU_noise_freq_dep.eps', 'epsc');

%% write out curve
out = [fne_vec', F_dot, sqrt(var_F_dot/N)];
T = array2table(out);
T.Properties.VariableNames(1:3) = {'f_ne','F_dot','err F_dot'}
writetable(T,['../data/sim_freq_dep_dg' num2str(dg) '_Dne' num2str(Dne) '_ts' num2str(ts)...
    '_N' num2str(N) '_K' num2str(K) '.csv']);

