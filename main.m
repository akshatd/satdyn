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
% x = [x, y, z, u, v, w, q_1, q_2, q_3, q_4, p, q, r]
x0 = zeros(13, 1);
x0(1) = planet.R + satellite.h; % total altitude
v_orbit = sqrt(planet.G * planet.M / norm(x0(1:3)));
x0(5) = v_orbit*cosd(satellite.i); % q
x0(6) = v_orbit*sind(satellite.i); % r
phi = 0; % rad angular positino about x
theta = 0; % rad angular position about y
psi = 0; % rad angular position about z
quat_scalar_first = eul2quat([psi, theta, phi], 'ZYX'); % quaternion
x0(7:10) = [quat_scalar_first(2:4), quat_scalar_first(1)]; % convert to scalar last
x0(11:13) = [0.01; 0.001; 0.02]; % angular velocity

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
surf(X, Y, Z, 'EdgeColor', 'none', 'DisplayName', 'Earth');
hold on;
axis equal;

% plot satellite
x_ode(:, 1:6) = x_ode(:, 1:6)/1000; % scale to km
plot3(x_ode(:, 1), x_ode(:, 2), x_ode(:, 3), 'r', 'LineWidth', 2, 'DisplayName', 'Satellite');
legend;

% plot angular position

quats = x_ode(:, 7:10);
euls = quat2eul([quats(:, 4), quats(:, 1:3)], 'ZYX'); % convert to scalar first then to euler angles

figure;
plot(t_ode, euls);
legend("$psi$", "$theta$", "$phi$", 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Angular Position');

function xdot = satellite_dynamics(t, x, planet, satellite)
% x = [x, y, z, u, v, w, q_1, q_2, q_3, q_4, p, q, r]

% translational dynamics
r_dot = x(4:6); % velocity
r = x(1:3); % position
r_norm = norm(r);
r_hat = r/r_norm;
Fg = -planet.G * planet.M * satellite.m / r_norm^2 * r_hat;
% add more forces here
F = Fg;
v_dot = F/satellite.m;

% rotational dynamics
q = x(7:10); % quaternion
omega = x(11:13); % angular velocity
q_dot = 0.5 * [
	q(4), -q(3), q(2);
	q(3), q(4), -q(1);
	-q(2), q(1), q(4);
	-q(1), -q(2), -q(3)
	] * omega;
% add more rotational dynamics here
omega_dot = [0; 0; 0];
xdot = [r_dot; v_dot; q_dot; omega_dot];
end