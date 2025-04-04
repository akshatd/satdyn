%% Satellite dynamics 6-DOF simulation
clc; clear; close all;
addpath("JpHelp\");

%% set up planet parameters
planet = {};
planet.M = 5.97219e24; % kg
planet.R = 6.63781e6; % m
planet.G = 6.67430e-11; % m^3/kg/s^2

%% Simulation parameters
Tend = 100; %[secs]
Ts = 3; % sampling period
Nsim = floor(Tend/Ts);

h0 = 500e3; % initial satellite height above ground [m]
i0 = 0; % initial satellite inclination [deg]

% initial conditions (see SatelliteClass for state and control definition)
% R0: initial position x,y,z in ECI frame 
R0 = [planet.R + h0;0;0]; % total altitude

% V0: initial velocity u,v,w in ECI frame
v_orbit = sqrt(planet.G * planet.M / norm(R0));
V0 = [0;v_orbit*cosd(i0);v_orbit*sind(i0)]; % pqr ??

% Q0: initial quaternion (Scaler-last) representing rotation of satellite body frame w.r.t. ECI frame
phi = 0*pi/180; % rad angular positino about x
theta = 0*pi/180; % rad angular position about y
psi = 0*pi/180; % rad angular position about z
quat_scalar_first = eul2quat([psi, theta, phi], 'ZYX'); % quaternion in inertial frame FA
Q0_A = [quat_scalar_first(2:4)'; quat_scalar_first(1)]; % convert to scalar last

% W0: initial angular velocity of satellite body w.r.t. ECI frame, expressed in body frame
W0 = [0.0; 0.0; 0.0];

%% Define satellite
sat_params = {};
sat_params.Nsim = Nsim;
sat_params.Ts = Ts;

sat_params.m = 9.755842;% Akshat 8; % mass (exclude reaction wheels) [kg]
sat_params.J_SatBody_C = [1.19*10^8 1.33*10^7 1.87*10^5;
                          1.33*10^7 8.79*10^7 -7.59*10^5;
                          1.87*10^5 -7.59*10^5 1.19*10^8]*10^-9; % Akshat diag([100;70;50])[kgm2] Inertia matrix of satellite body (exclude reaction wheels), wrt satellite's center of mass, resolved in body frame
sat_params.J_SatBody_C_Mean = mean(diag(sat_params.J_SatBody_C));

sat_params.I_ws = 0.5*0.137*(23/1000)^2 ;%Estimated from NanoAvio 0.5mr2 10;  % spin axis moments of inertia (wrt wheel's center-of-mass, in principal wheel frame) [kgm2]
% sat_params.gs_b_arr = [1 0 0;0 1 0; 0 0 1;1/sqrt(3) 1/sqrt(3) 1/sqrt(3);]'; % unit spin axis in satellite body's frame Fb [column vector x N_react] matrix
sat_params.gs_b_arr = [2 0 1; -2 0 1; 0 2 1; 0 -2 1]'/sqrt(5); % NanoAvio 4RW0 configuration
sat_params.gs_b_arr = [1 0 0;-1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]'; % 6 configuration
sat_params.MassReactionWheel = 0.137; %2; % mass of SINGLE reaction wheel, to be added to satellite mass [Kg]


% Initial reaction wheel speed
W_React0arr = zeros(size(sat_params.gs_b_arr,2),1);

% -1 = do nothing%
% 505 = Attitude rate control option 1
% 506 = 3D attitude control option 1 using PID
% 507 = 3D attitude control option 2 using basic nonlinear MPC
% 508 = 3D attitude control option 3 (P-outer, W-inner)
% 509 = 3D attitude control option 4 using reaction wheels nonlinear MPC


sat_params.ControlScheme = 506;  sat_params.Name = 'PID';
Satellite_PID = SatelliteClass(sat_params,R0,V0,Q0_A,W0,W_React0arr,planet);

sat_params.ControlScheme = 509;  sat_params.Name = 'NMPC1'; 
sat_params.q_err_quat_vec_BR_Body = 20;
sat_params.q_Delta_W_Body = 500;
sat_params.r_Delta_U = 0.00;
sat_params.r_U = 0.01;
sat_params.Renergy = 0;
Satellite_NMPC1 = SatelliteClass(sat_params,R0,V0,Q0_A,W0,W_React0arr,planet);

sat_params.ControlScheme = 509;  sat_params.Name = 'NMPC2'; 
sat_params.q_err_quat_vec_BR_Body = 20;
sat_params.q_Delta_W_Body = 1000;
sat_params.r_Delta_U = 0.00;
sat_params.r_U = 0.0075;
sat_params.Renergy = 10;
Satellite_NMPC2 = SatelliteClass(sat_params,R0,V0,Q0_A,W0,W_React0arr,planet);

SatelliteArr = [Satellite_PID,Satellite_NMPC1];%,Satellite_NMPC2];
nSat = length(SatelliteArr);
%% Simulate 
QuatRefA_A = Q0_A; % set initial quat refA

for k = 0:(Nsim-1)
    SimT = k*Ts;

    % set position commands 
    PosRef = zeros(3,1);
    VelRef = zeros(3,1);

    % Set attitude commands    
    if false % Cycle through various quaternions
        QuatRefA_A = zeros(4,1); QuatRefA_A(1 + floor(mod(k,200)/50)) = 1; %QuatRefA_Body(1) = 1;% 
    else % use a given profile
        GetQuatProfile;
    end
    % QuatRefA_Body = [zeros(3,1);1]; 

    % Set attitude rate reference
    if true % numerical differentiate QuatRefA_Body, followed by least-squares estimate of OmegaRefA_Body
        if k == 0 QuatRefA_A_Prev = QuatRefA_A;end
        delta_QuatRefA_A = QuatRefA_A - QuatRefA_A_Prev;
        Amat = Ts/2*[
	            QuatRefA_A(4), -QuatRefA_A(3), QuatRefA_A(2);
	            QuatRefA_A(3), QuatRefA_A(4), -QuatRefA_A(1);
	            -QuatRefA_A(2), QuatRefA_A(1), QuatRefA_A(4);
	            -QuatRefA_A(1), -QuatRefA_A(2), -QuatRefA_A(3)
	            ];
        OmegaRefA_R = (Amat'*Amat)\Amat'*delta_QuatRefA_A;
        O_RefA_K = quat2dcm([QuatRefA_A(end);QuatRefA_A(1:3)]');
        OmegaRefA_A = O_RefA_K'*OmegaRefA_R;
        QuatRefA_A_Prev = QuatRefA_A;
    else % Constant 
        OmegaRefA_A = zeros(3,1); %OmegaRefA_A(1) = 1*pi/180*sin(2*pi/1000*SimT); OmegaRefA_A(2) = 3*pi/180*cos(2*pi/1500*SimT); OmegaRefA_A(3) = 0.5*pi/180*sin(2*pi/750*SimT);
    end  

    % OmegaRefA_A = [0.00;0;-0.00]; QuatRefA_A = [zeros(3,1);1]; % do nothing

    for s = 1:nSat
        SatelliteArr(s).Step(PosRef,VelRef,QuatRefA_A,OmegaRefA_A);
    end
   
end

%% Plot results
close all
figN = 1;
for s = 1:nSat
    SatellitePlot = SatelliteArr(s); % select satellite to plot
    statesArr = SatellitePlot.statesArr(:,1:end-1);
    UarrStore = SatellitePlot.UarrStore;
    tArr = SatellitePlot.tArr;
    OmegaCmdA_BodyArr = SatellitePlot.OmegaCmdA_BodyArr;
    QuatRefA_A_Arr = SatellitePlot.QuatRefA_A_Arr;
    Err_BR_AngRadArr = SatellitePlot.Err_BR_AngRadArr;
    
    PlotAttitude1;

end


PlotAttitude2;
return;

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
xlabel('x'); ylabel('y'); zlabel('z');
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
legend("$\psi$", "$\theta$", "$\phi$", 'Interpreter', 'latex');
title('Angular Position');
ylabel('(deg)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,statesArr(11,:)*180/pi,'r',tArr,OmegaCmdA_BodyArr(1,:)*180/pi,'r--', 'LineWidth', 2);
plot(tArr,statesArr(12,:)*180/pi,'b',tArr,OmegaCmdA_BodyArr(2,:)*180/pi,'b--', 'LineWidth', 2);
plot(tArr,statesArr(13,:)*180/pi,'k',tArr,OmegaCmdA_BodyArr(3,:)*180/pi,'k--', 'LineWidth', 2);
legend("$p$","$p_{ref}$", "$q$", "$q_{ref}$", "$r$", "$r_{ref}$", 'Interpreter', 'latex');
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