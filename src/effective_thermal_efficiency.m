%EFFECTIVE_THERMAL_EFFICIENCY calculates the effective thermal efficiency
%and plots a comparison with Carnot efficiency
%
% OUTPUTS:
%  outputs eps figure of thermal efficiency compared with Carnot efficiency
%
% author:  JEhrich
% version: 1.3 (2022-12-12)
% changes: moved inset title to y label
clear
close all
clc
% set font size, line width, and marker size
fS = 21;
lW = 2.7;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% system parameters
Dne_vec = logspace(-2,10,1E3);
% sampling time
ts = 1/40;
% scaled effective mass
dg = 0.38;
% nonequilibrium noise frequency
fne = 1E4;


%% semi-analytic calculation of information power and free energy gain
% data structures
P_info = nan(length(Dne_vec),1);
F_dot = nan(length(Dne_vec),1);

for ii = 1:length(Dne_vec)
    Dne = Dne_vec(ii);
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
    P_info(ii) = (S_r_m - S_r)/ts;

    % output power
    F_dot(ii) = dg*sqrt(2)*sqrt(1 + Dne)*exp(-dg^2/(2 + 2*Dne))/(sqrt(pi)*(1 + erf(dg*sqrt(2)/(2*sqrt(1 + Dne)))));

end

% find Dne for which net output becomes positive
Dne_pos_output = Dne_vec(find(F_dot - P_info > 0,1));


%% plot comparison of thermal efficiency with Carnot efficiency


figure();
patch([min(Dne_vec),Dne_pos_output,Dne_pos_output,min(Dne_vec)],...
    [0,0,1,1],[1,1,1]*0.8,'EdgeColor','none');
hold on;
plot(Dne_vec,1-1./(1+Dne_vec),'k--','LineWidth',lW);
plot(Dne_vec,1-P_info./F_dot,'k-','LineWidth',lW);
set(gca,'XScale','log');
axis([min(Dne_vec),max(Dne_vec),0,1]);
set(gca,'FontSize',fS);
set(gca, 'Layer', 'top')
set(gca,"XTick",[1E-2, 1E0, 1E2, 1E4, 1E6, 1E8, 1E10])
xlabel('$D_\mathrm{ne}$','Interpreter','latex');
ylabel('efficiency','Interpreter','latex');
text(5E4,0.77,'$\eta_\mathrm{eff}$',...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fS,'Interpreter','latex');
text(2E0,0.9,'$\eta_\mathrm{C}$',...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fS,'Interpreter','latex');

% inset plotting difference
diff = 1-1./(1+Dne_vec') - (1-P_info./F_dot);
% thermal efficiency only meaningful for positive net output power
diff(F_dot-P_info < 0) = nan;
ax2 = axes('Position',[0.65,0.3,0.25,0.25]);
loglog(Dne_vec,diff,'k:','LineWidth',lW-0.3);
axis([1E2,max(Dne_vec),1E-4,1]);
set(gca,"XTick",[1E-2, 1E2, 1E6, 1E10])
set(gca,'FontSize',fS-4);
xlabel('$D_\mathrm{ne}$','Interpreter','latex');
ylabel('$\eta_\mathrm{C} - \eta_\mathrm{eff}$','Interpreter','latex');

% save
saveas(gcf, '../doc/thermal_efficiency.eps','epsc')

