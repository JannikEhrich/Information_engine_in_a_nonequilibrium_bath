%EFFICIENCY_COMPARISON calculates information power required to run
%information engine in a nonequilibrium bath and plots comparison of
%efficiency with conventional engine
%
% OUTPUTS:
%  outputs eps figure of information power and efficiency
%
% author:  JEhrich
% version: 1.2 (2022-07-06)
% changes: chaged ts to 1/42 to align with experiments
clear
close all
clc
% set font size, line width, and marker size
fS = 18;
lW = 2.7;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% system parameters
Dne_vec = logspace(-2,4,20);
% sampling time
ts = 1/42;
% scaled effective mass
dg = 0.38;
% nonequilibrium noise frequency
fne = 1E4;

% maximum Dne in experiment
Dne_exp_max = 1E2;

%% simulation parameters
% number of samples
K = 1E6;
% number of steps to equilibrate
Kini = 1E4;
% number of bins for histogram
n_bins = 1E2;

%% simulate information ratchet
% data structures
v = nan(length(Dne_vec),1);
F_dot = nan(length(Dne_vec),1);
P_trap = nan(length(Dne_vec),1);
P_info = nan(length(Dne_vec),1);

tic
for ii = 1:length(Dne_vec)
    ii
    Dne = Dne_vec(ii);
    % reach steady-state ratcheting
    [x_traj, l_traj, zeta_traj] = sim_OU_ratchet(dg,Dne,fne,ts,Kini);
    l0 = l_traj(end);
    % simulate a steady-state trajectory
    [x_traj, l_traj, ~] = sim_OU_ratchet(dg,Dne,fne,ts,K,...
        x_traj(end),l0,zeta_traj(end));
    % velocity
    v(ii) = (l_traj(end) - l0)/K/ts;
    % rate of free energy gain
    F_dot(ii) = dg*(l_traj(end) - l0)/K/ts;
    % trap power
    P_trap(ii) = mean(1/2*(x_traj(2:end)-l_traj(2:end)).^2 ...
        - 1/2*(x_traj(2:end)-l_traj(1:end-1)).^2);

    % relative distribution at same time and entropy
    [c_same, bins_s] = hist(l_traj - x_traj,n_bins);
    ker = log(c_same/K/diff(bins_s(1:2)));
    ker(isnan(ker) | isinf(ker)) = 0;
    S_same = -sum(c_same/K.*ker);
    % relative distribution at differnt times and entropy
    [c_diff, bins_d] = hist(l_traj(1:end-1) - x_traj(2:end),n_bins);
    ker = log(c_diff/(K-1)/diff(bins_d(1:2)));
    ker(isnan(ker) | isinf(ker)) = 0;
    S_diff = -sum(c_diff/(K-1).*ker);
    % minimum information power
    P_info(ii) = (S_diff - S_same)/ts;

    % save initial and final histograms
    if ii == 1
        pi_r_ini = c_same/K/diff(bins_s(1:2));
        r_ini = bins_s;
        pi_rm_ini = c_diff/(K-1)/diff(bins_d(1:2));
        rm_ini = bins_d;
    end
    if ii == 13
        Dne_mid = Dne;
        pi_r_mid = c_same/K/diff(bins_s(1:2));
        r_mid = bins_s;
        pi_rm_mid = c_diff/(K-1)/diff(bins_d(1:2));
        rm_mid = bins_d;
    end
    if ii == length(Dne_vec)
        pi_r_fin = c_same/K/diff(bins_s(1:2));
        r_fin = bins_s;
        pi_rm_fin = c_diff/(K-1)/diff(bins_d(1:2));
        rm_fin = bins_d;
    end
end
toc

%% semi-analytic calculation of information power
% vector for analytical evaluation
Dne_ana = logspace(log10(min(Dne_vec)),log10(max(Dne_vec)),1E3);

% data structure
P_info_ana = nan(length(Dne_ana),1);

for ii = 1:length(Dne_ana)
    Dne = Dne_ana(ii);
    % grids for stationary distributions
    r = linspace(0,4*sqrt(1+Dne),1E4);
    dr = diff(r(1:2));
    r_m = linspace(-1*sqrt(1+Dne),4*sqrt(1+Dne),1E4);
    dr_m = diff(r_m(1:2));
    % normalization
    A = sqrt(2/pi/(1+Dne))/(1+ erf(dg/sqrt(2*(1+Dne))));
    % stationary distributions
    pi_r = A*exp(-(r-dg).^2/2/(1+Dne));
    pi_r_m = -A*exp(-(dg - r_m).^2/(2 + 2*Dne))...
        .*(erf(sqrt(2)*((dg - r_m)*exp(-ts) - dg)/(2*sqrt(1 - exp(-2*ts))*sqrt(1 + Dne))) - 1)/2;
    % entropies
    ker = log(pi_r);
    ker(isnan(ker) | isinf(ker)) = 0;
    S_r = -sum(pi_r.*ker)*dr;
    ker = log(pi_r_m);
    ker(isnan(ker) | isinf(ker)) = 0;
    S_r_m = -sum(pi_r_m.*ker)*dr_m;
    % information power
    P_info_ana(ii) = (S_r_m - S_r)/ts;

    % output power
    F_dot(ii) = dg*sqrt(2)*sqrt(1 + Dne)*exp(-dg^2/(2 + 2*Dne))/(sqrt(pi)*(1 + erf(dg*sqrt(2)/(2*sqrt(1 + Dne)))));


    % save initial, middle, and final distributions
    if ii == 1
        pi_r_ini_ana = pi_r;
        r_ini_ana = r;
        pi_rm_ini_ana = pi_r_m;
        rm_ini_ana = r_m;
    end
    % find middle index
    [~,ind_mid] = min(abs(Dne_mid-Dne_ana)); 
    if ii == ind_mid
        pi_r_mid_ana = pi_r;
        r_mid_ana = r;
        pi_rm_mid_ana = pi_r_m;
        rm_mid_ana = r_m;
    end
    if ii == length(Dne_ana)
        pi_r_fin_ana = pi_r;
        r_fin_ana = r;
        pi_rm_fin_ana = pi_r_m;
        rm_fin_ana = r_m;
    end
end



%% plot information power 
figure('Position',[400,1000,560,490]);
ax1 = axes('Position',[0.13 0.51 0.82 0.46]);
semilogx(Dne_ana,P_info_ana,'k-','LineWidth',lW,'MarkerSize',mS);
hold on;
semilogx(Dne_vec,P_info,'ks','LineWidth',lW,'MarkerSize',mS);
xlabel('$D_\mathrm{ne}$','Interpreter','latex');
ylabel('$P_\mathrm{info}$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'approximation','numerics'},...
    'Location','northwest');
set(gca,'YLim',[4.5,6.7],'XLim',[min(Dne_vec),max(Dne_vec)]);
text(min(Dne_vec)/7.5,6.7,'(a)','Interpreter','latex','FontSize',fS);

% plot distibutions
ax2 = axes('Position',[0.13 0.09 0.22 0.26]);
plot(nan,nan,'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(nan,nan,'b','LineWidth',lW,'MarkerSize',mS);
bar(r_ini,pi_r_ini,'r','EdgeColor','none','BarWidth',1,'FaceAlpha',0.3);
box off;
bar(rm_ini,pi_rm_ini,'b','EdgeColor','none','BarWidth',1,'FaceAlpha',0.3);
plot(rm_ini_ana,pi_rm_ini_ana,'b','LineWidth',lW,'MarkerSize',mS);
plot(r_ini_ana,pi_r_ini_ana,'r','LineWidth',lW,'MarkerSize',mS);
set(gca,'XLim',[min(rm_ini_ana),max(rm_ini_ana)]);
set(gca,'Ylim',[0,0.7])
set(gca,'FontSize',fS-2);
xlabel('$r$','Interpreter','latex');
%title(['$D_\mathrm{ne}=10^{' num2str(log10(Dne_vec(1))) '}$'],'FontWeight','normal','Interpreter','latex','FontSize',fS-2)
annotation('textbox',[.17 .18 .1 .2], ...
    'String',['$D_\mathrm{ne}=10^{' num2str(log10(Dne_vec(1))) '}$'],...
    'EdgeColor','none','Interpreter','latex','FontSize',fS-2)
text(-3.7,0.7,'(b)','Interpreter','latex','FontSize',fS);


ax3 = axes('Position',[0.43 0.09 0.22 0.26]);
plot(nan,nan,'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(nan,nan,'b','LineWidth',lW,'MarkerSize',mS);
bar(r_mid,pi_r_mid,'r','EdgeColor','none','BarWidth',1,'FaceAlpha',0.3);
box off;
bar(rm_mid,pi_rm_mid,'b','EdgeColor','none','BarWidth',1,'FaceAlpha',0.3);
plot(rm_mid_ana,pi_rm_mid_ana,'b','LineWidth',lW,'MarkerSize',mS);
plot(r_mid_ana,pi_r_mid_ana,'r','LineWidth',lW,'MarkerSize',mS);
set(gca,'XLim',[min(rm_mid_ana),max(rm_mid_ana)]);
set(gca,'Ylim',[0,0.11])
set(gca,'FontSize',fS-2);
xlabel('$r$','Interpreter','latex');
annotation('textbox',[.47 .18 .1 .2], ...
    'String',['$D_\mathrm{ne}=' num2str(round(Dne_mid)) '$'],...
    'EdgeColor','none','Interpreter','latex','FontSize',fS-2)

ax4 = axes('Position',[0.73 0.09 0.22 0.26]);
plot(nan,nan,'r','LineWidth',lW,'MarkerSize',mS);
hold on;
plot(nan,nan,'b','LineWidth',lW,'MarkerSize',mS);
h = bar(r_fin,pi_r_fin,'r','EdgeColor','none','BarWidth',1,'FaceAlpha',0.3);
% remove scientific notation from y-axis
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
xtickformat('%.0f')
box off;
bar(rm_fin,pi_rm_fin,'b','EdgeColor','none','BarWidth',1,'FaceAlpha',0.3);
plot(rm_fin_ana,pi_rm_fin_ana,'b','LineWidth',lW,'MarkerSize',mS);
plot(r_fin_ana,pi_r_fin_ana,'r','LineWidth',lW,'MarkerSize',mS);
set(gca,'XLim',[min(rm_fin_ana),max(rm_fin_ana)]);
set(gca,'Ylim',[0,9E-3])
set(gca,'FontSize',fS-2);
xlabel('$r$','Interpreter','latex');
annotation('textbox',[.77 .18 .1 .2], ...
    'String',['$D_\mathrm{ne}\!=\!10^{' num2str(log10(Dne_vec(end))) '}$'],...
    'EdgeColor','none','Interpreter','latex','FontSize',fS-2)

% save
saveas(gcf, '../doc/information_power.eps','epsc')

%% plot output, input power and efficiency
% velocity from output power
v = F_dot/dg;

% input and output power
figure('Position',[400,1000,560,600]);
ax1 = axes('Position',[0.13 0.55 0.77 0.42]);
loglog(Dne_ana(Dne_ana<Dne_exp_max),F_dot(Dne_ana<Dne_exp_max),'LineWidth',20,'MarkerSize',mS,'color',0.8*[1,1,1]);
hold on;
loglog(Dne_ana(Dne_ana<Dne_exp_max),P_info_ana(Dne_ana<Dne_exp_max),'LineWidth',20,'MarkerSize',mS,'color',[1,0.8,0.8]);
loglog(Dne_ana,F_dot,'k','LineWidth',lW,'MarkerSize',mS);
loglog(Dne_ana,P_info_ana,'r','LineWidth',lW,'MarkerSize',mS);
set(gca,'XTickLabel',[]);
ylabel('power','Interpreter','latex');
set(gca,'FontSize',fS);
set(gca,'XLim',[min(Dne_vec),max(Dne_vec)]);
set(gca,'YLim',[1E-1, 3.5E1]);
text(1E1,0.52,'$\dot F$','Interpreter','latex','FontSize',fS);
text(1E-1,7.5E0,'$P_\mathrm{info}$','Interpreter','latex','FontSize',fS,'color','r');
text(1E-3,3.5E1,'\textbf{(a)}','Interpreter','latex','FontSize',fS);

% efficiency
ax2 = axes('Position',[0.13 0.085 0.77 0.42]);
semilogx(Dne_ana(Dne_ana<Dne_exp_max),F_dot(Dne_ana<Dne_exp_max)./P_info_ana(Dne_ana<Dne_exp_max),'LineWidth',20,'MarkerSize',mS,'color',0.8*[1,1,1]);
hold on;
semilogx(Dne_ana,F_dot./P_info_ana,'k','LineWidth',lW,'MarkerSize',mS);
xlabel('$D_\mathrm{ne}$','Interpreter','latex');
ylabel('efficiency','Interpreter','latex');
set(gca,'FontSize',fS);
set(gca,'XLim',[min(Dne_vec),max(Dne_vec)]);
set(gca,'YLim',[-0.04,1.5]);
text(3E2,0.8,'$\dot F/P_\mathrm{info}$','Interpreter','latex','FontSize',fS,'color','k');
text(1E-3,1.5,'\textbf{(b)}','Interpreter','latex','FontSize',fS);

% inset
ax3 = axes('Position',[0.22 0.24 0.33 0.23]);
semilogx(F_dot,ones(size(F_dot)),'k--','LineWidth',lW-0.5);
hold on;
loglog(F_dot(Dne_ana<Dne_exp_max),F_dot(Dne_ana<Dne_exp_max)./P_info_ana(Dne_ana<Dne_exp_max),'LineWidth',20,'MarkerSize',mS,'color',[1,0.8,0.8]);
semilogx(F_dot,F_dot./P_info_ana,'r-','LineWidth',lW-0.5);
semilogx(F_dot,F_dot./(v.^2 + v*dg),'b-','LineWidth',lW-0.5);
set(gca,'XLim',[min(F_dot),max(F_dot)]);
set(gca,'YLim',[-0.08,1.5]);
text(1E0,1.2,'(ratchet)','Interpreter','latex','FontSize',fS-2,'color','r');
text(5E0,0.2,'(pulling)','Interpreter','latex','FontSize',fS-2,'color','b');
set(gca,'FontSize',fS-2);
set(gca,'XTick',[1E-1,1E0,1E1]);
xlabel('$\dot F$','Interpreter','latex')

saveas(gcf, '../doc/compare_info_efficiency.eps','epsc')

% write csv of information power data
out = [Dne_ana',F_dot,P_info_ana];
T = array2table(out);
T.Properties.VariableNames(1:3) = {'Dne','F_dot','P_info'};
writetable(T,['../data/info_power_dg_' num2str(dg) '_ts_' num2str(ts)...
    '_fne_' num2str(fne) '.csv']);

% write csv of efficiency inset data
out = [F_dot, F_dot./P_info_ana, F_dot./(v.^2 + v*dg)];
T = array2table(out);
T.Properties.VariableNames(1:3) = {'F_dot','F_dot/P_info','F_dot/P_trap (conventional)'};
writetable(T,['../data/info_efficiency_inset_dg_' num2str(dg) '_ts_' num2str(ts)...
    '_fne_' num2str(fne) '.csv']);
