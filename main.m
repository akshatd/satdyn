%% Satellite dynamics 6-DOF simulation
clc; clear; close all;

% set up planet and satellite parameters
planet = {};
planet.M = 5.97219e24; % kg
planet.R = 6.63781e6; % m
planet.G = 6.67430e-11; % m^3/kg/s^2

satellite = {};
satellite.m = 8; % kg
satellite.h = 500e3; % m
satellite.i = 0; % deg, inclination

% initial conditions
x0 = zeros(6, 1);
x0(1) = planet.R + satellite.h; % total altitude
v_orbit = sqrt(planet.G * planet.M / norm(x0(1:3)));
x0(5) = v_orbit*cosd(satellite.i); % q
x0(6) = v_orbit*sind(satellite.i); % r

% simulation
tspan = [0, 10000]; % s
[t_ode, x_ode] = ode45(@(t, x) satellite_dynamics(t, x, planet, satellite), tspan, x0);

% plot earth
[X, Y, Z] = sphere();
% scale to earth radius in km
X = X * planet.R/1000;
Y = Y * planet.R/1000;
Z = Z * planet.R/1000;
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
hold on;
axis equal;
% view(3);

% plot satellite
x_ode = x_ode/1000; % scale to km
plot3(x_ode(:, 1), x_ode(:, 2), x_ode(:, 3), 'r', 'LineWidth', 2, 'DisplayName', 'Satellite');
legend;



function xdot = satellite_dynamics(t, x, planet, satellite)
% x = [x, y, z, vx, vy, vz]
r_dot = x(4:6); % velocity
r = x(1:3); % position
r_norm = norm(r);
r_hat = r/r_norm;
Fg = -planet.G * planet.M * satellite.m / r_norm^2 * r_hat;
% add more forces here
F = Fg;
v_dot = F/satellite.m;
xdot = [r_dot; v_dot];
end