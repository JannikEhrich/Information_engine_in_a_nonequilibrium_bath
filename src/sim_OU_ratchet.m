function [x_traj,l_traj,zeta_traj] = sim_OU_ratchet(dg,Dne,fne,ts,K,x0,l0,zeta0)
%SIM_OU_ratchet simulates info ratchet dynamics with additional
%Ornstein-Uhlenbeck noise with correlation frequency fne and strength Dne
%
% INPUTS:
%      dg: scaled effective mass
%     Dne: active O-U noise stength
%     fne: active O-U noise correlation frequency
%      ts: sampling time
%       K: number of samples
%      x0: (optional) initial particle position
%      l0: (optional) initial ratchet position
%   zeta0: (optional) initial noise value
%
% OUTPUTS:
%      x_traj: particle trajectory
%      l_traj: trap trajectory
%   zeta_traj: noise trajectory
%
% author:  JEhrich
% version: 1.4 (2022-06-30)
% changes: removed variable feedbakc gain                 

%% initialize positions with equilibrium values
if nargin < 8
    x0 = rand*sqrt(1+Dne)-dg;
    l0 = 0;
    zeta0 = rand*sqrt(Dne*fne);
end
x = x0;
l = l0;
zeta = zeta0;

%% initialize trajectory vectors
x_traj = nan(K+1,1);
l_traj = nan(K+1,1);
zeta_traj = nan(K+1,1);
x_traj(1) = x;
l_traj(1) = l;
zeta_traj(1) = zeta;

%% variances for stochastic dynamics
c_xx = (4*Dne*fne^2*exp(-(1 + fne)*ts) - Dne*fne*(1 + fne)*exp(-2*fne*ts) - (1 + (1 + Dne)*fne^2 - 2*fne)*(1 + fne)*exp(-2*ts) + (1 + (1 + Dne)*fne)*(fne - 1)^2)/((fne - 1)^2*(1 + fne));
c_xzeta = fne*Dne*(-2*fne*exp(-(1 + fne)*ts) + (1 + fne)*exp(-2*fne*ts) + fne - 1)/(fne^2 - 1);
c_zetazeta = -Dne*fne*(exp(-2*fne*ts) - 1);

%% main loop
for ii = 1:K
    % new particle position and active noise
    mu_x = (-zeta*exp(-fne*ts) + ((dg - l + x)*fne - dg + l - x + zeta)*exp(-ts) + (fne - 1)*(l - dg))/(fne - 1);
    mu_zeta = zeta*exp(-fne*ts);
    % draw next zeta state
    zeta = mu_zeta + sqrt(c_zetazeta)*randn;
    % draw next x state
    x = mu_x + c_xzeta/c_zetazeta*(zeta - mu_zeta)...
        + sqrt(c_xx - c_xzeta^2/c_zetazeta)*randn;
    % feedback
    if x > l
        l = l + 2*(x-l);
    end
    % write trajectory vectors
    x_traj(ii+1) = x;
    l_traj(ii+1) = l;
    zeta_traj(ii+1) = zeta;
end


end