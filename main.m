%% Satellite dynamics 6-DOF simulation
clc; clear; close all;

%% set up planet parameters
planet = {};
planet.M = 5.97219e24; % kg
planet.R = 6.63781e6; % m
planet.G = 6.67430e-11; % m^3/kg/s^2

%% Simulation parameters
Tend = 10000; %[secs]
Ts = 10; % sampling period
Nsim = Tend/Ts;

h0 = 500e3; % initial satellite height above ground [m]
i0 = 0; % initial satellite inclination [deg]

% initial conditions (see SatelliteClass for state and control definition)
% R0: initial position x,y,z in ECI frame 
R0 = [planet.R + h0;0;0]; % total altitude

% V0: initial velocity u,v,w in ECI frame
v_orbit = sqrt(planet.G * planet.M / norm(R0));
V0 = [0;v_orbit*cosd(i0);v_orbit*sind(i0)]; % pqr ??

% Q0: initial quaternion (Scaler-last) representing rotation of satellite body frame w.r.t. ECI frame
phi = 0; % rad angular positino about x
theta = 0; % rad angular position about y
psi = 0; % rad angular position about z
quat_scalar_first = eul2quat([psi, theta, phi], 'ZYX'); % quaternion
Q0 = [quat_scalar_first(2:4)'; quat_scalar_first(1)]; % convert to scalar last

% W0: initial angular velocity of satellite body w.r.t. ECI frame, expressed in body frame
W0 = [0.0; 0.0; 0.0];

%% Define satellite
sat_params = {};
sat_params.Nsim = Nsim;
sat_params.Ts = Ts;

sat_params.m = 8; % mass (exclude reaction wheels) [kg]
sat_params.J_SatBody_C = diag([100;70;50]); % TO SET! [kgm2] Inertia matrix of satellite body (exclude reaction wheels), wrt satellite's center of mass, resolved in body frame

sat_params.I_ws = 10;  % TO SET! spin axis moments of inertia (wrt wheel's center-of-mass, in principal wheel frame) [kgm2]
sat_params.gs_b_arr = [1 0 0;0 1 0; 0 0 1;1/sqrt(3) 1/sqrt(3) 1/sqrt(3);]'; % unit spin axis in satellite body's frame Fb [column vector x N_react] matrix
sat_params.MassReactionWheel = 2; % mass of SINGLE reaction wheel, to be added to satellite mass [Kg]

sat_params.ControlScheme = -1; % do nothing

% Initial reaction wheel speed
W_React0arr = zeros(size(sat_params.gs_b_arr,2),1);

sat_params.Name = 'Satellite1';
Satellite1 = SatelliteClass(sat_params,R0,V0,Q0,W0,W_React0arr,planet);

sat_params.Name = 'Satellite2';
Satellite2 = SatelliteClass(sat_params,R0,V0,Q0,W0,W_React0arr,planet);
%% Simulate 

for k = 0:(Nsim-1)
    % set commands 
    PosRef = zeros(3,1);
    VelRef = zeros(3,1);
    QuatRef = zeros(3,1);
    OmegaRef = zeros(3,1);

    Satellite1.Step(PosRef,VelRef,QuatRef,OmegaRef);

    Satellite2.Step(PosRef,VelRef,QuatRef,OmegaRef); % currently redundant as it is identical to satellite 1
end

%% Plot results
figN = 1;
SatellitePlot = Satellite1; % select satellite to plot
statesArr = SatellitePlot.statesArr(:,1:end-1);
UarrStore = SatellitePlot.UarrStore;
tArr = SatellitePlot.tArr;

% plot earth and satellite in 3D
[X, Y, Z] = sphere();
% scale to earth radius in km
X = X * planet.R/1000;
Y = Y * planet.R/1000;
Z = Z * planet.R/1000;
figure(figN); figN = figN + 1;
surf(X, Y, Z, 'EdgeColor', 'none', 'DisplayName', 'Earth');
hold on;
axis equal;
plot3(statesArr(1,:)/1000, statesArr(2,:)/1000, statesArr(3,:)/1000, 'r', 'LineWidth', 2, 'DisplayName', 'Satellite');
legend;

% plot satellite states
figure(figN); figN = figN + 1; lm = []; kk = 1;
nRows = 3; nCols = 2;

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,statesArr(1,:)/1000,tArr,statesArr(2,:)/1000,tArr,statesArr(3,:)/1000, 'LineWidth', 2)
legend('x','y','z')
title('Position')
ylabel('(km)')

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,statesArr(4,:)/1000,tArr,statesArr(5,:)/1000,tArr,statesArr(6,:)/1000, 'LineWidth', 2)
legend('u','v','w')
title('velocity')
ylabel('(km/s)')

quats = statesArr(7:10,:)';
euls = quat2eul([quats(:, 4), quats(:, 1:3)], 'ZYX'); % convert to scalar first then to euler angles

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,euls*180/pi, 'LineWidth', 2)
legend("$psi$", "$theta$", "$phi$", 'Interpreter', 'latex');
title('Angular Position');
ylabel('(deg)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,statesArr(11,:)*180/pi,tArr,statesArr(12,:)*180/pi,tArr,statesArr(13,:)*180/pi, 'LineWidth', 2)
legend("$p$", "$q$", "$r$", 'Interpreter', 'latex');
title('Angular velocity');
ylabel('(deg/s)');

% Plot reaction wheels
lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
leg = [];
title('Reaction wheel angular velocity');
for r = 1:SatellitePlot.params.N_react
    leg(r) = plot(tArr,statesArr(13+r,:)*180/pi, 'LineWidth', 2,'DisplayName',num2str(r));   
end
ylabel('(deg/s)');
xlabel('Time (s)');
legend(leg, 'Interpreter', 'latex');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
leg = [];
title('Reaction wheel angular acceleration');
for r = 1:SatellitePlot.params.N_react
    leg(r) = plot(tArr,UarrStore(r,:)*180/pi, 'LineWidth', 2,'DisplayName',num2str(r));   
end
ylabel('(deg/s2)');
xlabel('Time (s)');
legend(leg, 'Interpreter', 'latex');

linkaxes(lm,'x');