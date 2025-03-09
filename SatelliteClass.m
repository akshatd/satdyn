classdef SatelliteClass < handle & matlab.mixin.Heterogeneous
    % Implements the satellite object
    % State definition: x = [x, y, z, u, v, w, q_1, q_2, q_3, q_4, p, q, r,OmeReac_1; OmeReac_2; ... ;OmeReac_N_react]
    % Controls definition: U = [OmeReacDot_1; OmeReacDot_2; ...; OmeReacDot_N_react];

    properties
        % User-provide
        params
        planet
        
        
        % derived
        statesArr % see above state definition
        SimCnt = 0; % simulation step count, starts from 0
        SimT = 0;   % simulation time = SimCnt*Ts
        UarrStore 
        tArr

    end
    
    methods (Access = public)
        function obj = SatelliteClass(params,R0,V0,Q0,W0,W_React0arr, planet)
            % inputs (vectors are columns):

            % params shall contain the following fields:
            % m = mass (exclude reaction wheels) [kg]
            % Nsim = number of simulation steps
            % Ts = sampling period
            % J_SatBody_C =  [kgm2] Inertia matrix of satellite body (exclude reaction wheels), wrt satellite's center of mass, resolved in body frame
             % N_react x reaction wheel (all are identical, except for wheel's orientation in satellite)
            % gs_b_arr    % unit spin axis in satellite body's frame Fb [column vector x N_react] matrix
            % I_ws    % spin axis moments of inertia (wrt wheel's center-of-mass, in principal wheel frame) [kgm2]
            % MassReactionWheel   % mass of SINGLE reaction wheel, to be added to satellite mass [Kg]
            % ControlScheme % integer

            % R0: initial position x,y,z in ECI frame 
            % V0: initial velocity u,v,w in ECI frame
            % Q0: initial quaternion (Scaler-last) representing rotation of satellite
            % body frame w.r.t. ECI frame
            % W0: initial angular velocity of satellite body w.r.t. ECI
            % frame, expressed in body frame
            % W_React0arr: initial angular velocity of reaction wheels
            
           
            % -------- derived params -------- 
            % N_react     % number of reaction wheel
            % gg_b_arr    % points in the wheel gimbal direction (just need to be perpendicular to gs_b for reaction wheels)
            % gt_b_arr    % completes the gimbal frame
            % Itot        % Total moments of inertia = satellite + reaction wheels


            % init states storage for entire simulation
          

            obj.params = params;
            obj.planet = planet;

            % Define initial state vector
            states0 = [R0;V0;Q0;W0;W_React0arr];

            obj.statesArr = zeros(length(states0),params.Nsim + 1);
            obj.statesArr(:,1) = states0;
            obj.tArr = zeros(1,params.Nsim);

            % Reaction wheel
            obj.params.N_react = size(obj.params.gs_b_arr,2);
            fprintf('%d reaction wheels \n',obj.params.N_react);

            for r = 1:obj.params.N_react 
                % Set the irrelevant gg_b_arr and gt_b_arr
                Z = null(obj.params.gs_b_arr(:,r)');
                obj.params.gg_b_arr(:,r) = Z(:,1); % Arbitrary assignment of gt and gg
                obj.params.gt_b_arr(:,r) = Z(:,2); 
            end

            obj.params.Itot = obj.params.J_SatBody_C; % TODO!! This is wrong
        end

        function U = GetControls(obj,PosRef,VelRef,QuatRef,OmegaRef)
            switch obj.params.ControlScheme
                otherwise
                    U = zeros(obj.params.N_react,1);
            end
        end

        function Step(obj,PosRef,VelRef,QuatRef,OmegaRef)
            % 1 step in the current simulation

            % Compute controls
            U = obj.GetControls(PosRef,VelRef,QuatRef,OmegaRef);

            if obj.SimCnt == 0 % init UarrStore based on number of controls
                obj.UarrStore = zeros(length(U),obj.params.Nsim);
            end

            % propagate state forward 1 step
            odeFun = @(t,x) obj.satellite_dynamics(t, x,U, obj.planet, obj.params);
            tspan = [obj.SimT (obj.SimT+obj.params.Ts)]; 
            xCurrent = obj.statesArr(:,obj.SimCnt + 1);
            [T,X] = ode45(odeFun, tspan, xCurrent); 

           % Store memory and propagate counter
           obj.UarrStore(:,obj.SimCnt + 1) = U;
           obj.statesArr(:,obj.SimCnt + 2) = X(end,:)';
           obj.tArr(:,obj.SimCnt + 1) = obj.SimT;
           obj.SimCnt = obj.SimCnt + 1;
           obj.SimT = obj.SimCnt*obj.params.Ts;

        end
        
       
    end % public

    methods (Static)
        function xdot = satellite_dynamics(t, x,U, planet, sat_params)
            % Compute the state's continous time-derivative
            % See claas definition for X and U definition
            
            % translational dynamics
            massTotal = sat_params.m + sat_params.N_react*sat_params.MassReactionWheel;
            r_dot = x(4:6); % velocity
            r = x(1:3); % position
            r_norm = norm(r);
            r_hat = r/r_norm;
            Fg = -planet.G * planet.M * massTotal / r_norm^2 * r_hat;
            % add more forces here
            F = Fg;
            v_dot = F/massTotal;
            
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

            % Reaction wheel dynamics
            omega_reac = x(14:(13+sat_params.N_react));
            external_torque = zeros(3,1); % Edit here for external torques
            I_wdot = -cross(omega,sat_params.Itot*omega) + external_torque;
            OmeReacDot = U(1:sat_params.N_react);
            for r = 1:sat_params.N_react % Sum impact of each reaction wheel
                % Compute omega_G = [ws;wt;wg] = angular velocity of
                % satellite relative to inertial, resolved in gimbal frame

                % get O_BG
                O_BG = [sat_params.gs_b_arr(:,r) sat_params.gt_b_arr(:,r) sat_params.gg_b_arr(:,r)];
                omega_G = O_BG'*omega; 

                I_wdot = I_wdot - sat_params.I_ws*OmeReacDot(r)*sat_params.gs_b_arr(:,r) ...
                                - sat_params.I_ws*omega_reac(r)*omega_G(3)*sat_params.gt_b_arr(:,r) ...
                                + sat_params.I_ws*omega_reac(r)*omega_G(2)*sat_params.gg_b_arr(:,r);

            end
            omega_dot = sat_params.Itot\I_wdot;
            xdot = [r_dot; v_dot; q_dot; omega_dot;OmeReacDot];
        end
    end % static
end

