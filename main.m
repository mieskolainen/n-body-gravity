% N-body gravity simulation with adaptive multilevel time stepping
%
% Run by [F5]
%
%
% Computational Astrophysics Course Project, ETH Zurich, 2012
% mikael.mieskolainen@cern.ch

close all;
clear;

%f = '2body';    % 2-body
%f = 'cemb';     % Circular equal mass binary
%f = 'eqmb';     % Eccentric equal mass binary
%f = 'cnemb';    % Circular non-equal mass binary
%f = 'enmb';     % Eccentric non-equal mass binary
%f = 'ko';       % Keplerian orbit
%f = '3bodyem';  % 3-body equal mass
f = '3bodynem';  % 3-body non-equal mass
%f = '4bodyem';  % 4-body equal mass

inputpath = 'input/';
filename = sprintf("%s%s", inputpath, f);

tmax = 3.0;      % Runtime
eps = 0.0001;    % Softening
n = 4;           % Number of timestep levels/bins


% Integrator method (different methods for different problems)
%method = 'Euler';
%method = 'DKD';
method = 'KDK';
%method = 'RK4';

% N-body integration
tic
[id,mass,pos,vel,dE_E,E,t,timesteps] = nbody(filename, tmax, n, eps, method);
toc


% Plots

bodies = length(id);

figure; subplot(2,2,1);
% Plot all particle trajectories
colors = {'b','r','g','k', ...
          'b','r','g','k', ...
          'b','r','g','k', ...
          'b','r','g','k'};

for i = 1:bodies
    plot3(squeeze(pos(i,1,:)), ...
          squeeze(pos(i,2,:)), squeeze(pos(i,3,:)), colors{i}); hold on;
end
view(3); axis equal; axis equal;
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');
title('Trajectories in $R^3$','interpreter','latex');

subplot(2,2,2);
for i = 1:bodies
    plot(squeeze(pos(i,1,:)), squeeze(pos(i,2,:)), colors{i}); hold on;
end
view(2); axis equal;
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
title('Trajectories in $R^2$','interpreter','latex');


subplot(2,2,3);
plot(t,dE_E);
xlabel('$t$ (sec)','interpreter', 'latex');
ylabel('$\Delta E_{tot}$','interpreter','latex');

subplot(2,2,4);
plot(t,E);
xlabel('$t$ (sec)','interpreter', 'latex');
ylabel('$E$','interpreter','latex');
l = legend('E_{tot}', 'E_k', 'E_p'); set(l,'fontsize',6,'location','southeast');


figure;
plot(t(1:end-1),timesteps');
xlabel('$t$ (sec)','interpreter', 'latex');
ylabel('$\Delta t_i$', 'interpreter', 'latex');
title('Timesteps of particles (color = particle)','interpreter','latex');

