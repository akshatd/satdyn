classdef SatelliteClass < handle & matlab.mixin.Heterogeneous
    % Implements the satellite object
    % State definition: x = [x, y, z, u, v, w, q_1, q_2, q_3, q_4, p, q, r,OmeReac_1; OmeReac_2; ... ;OmeReac_N_react]
    % note that [q1 q2 q3 q4] = quat_BA_(B or A), the quaternion
    % representing rotation of body FB from inertial FA, resolved in FB or
    % FA (same result)
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

        OmegaRefA_BodyArr   % store OmegaRefA_Body
        QuatRefA_BodyArr   % store QuatRefA_Body

        OmegaRef_Prev = zeros(3,1);

        Body_int_506 = zeros(3,1); % integrator vector state for option 506
        Delta_W0_Body % store initial delta W

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
            % Itot_Body        % Total moments of inertia = satellite + reaction
            % wheels, of full body, wrt full body's center of mass,
            % expressed in body frame


            % init states storage for entire simulation
          

            obj.params = params;
            obj.planet = planet;

            % Define initial state vector
            states0 = [R0;V0;Q0;W0;W_React0arr];

            obj.statesArr = zeros(length(states0),params.Nsim + 1);
            obj.statesArr(:,1) = states0;
            obj.tArr = zeros(1,params.Nsim);
            obj.OmegaRefA_BodyArr = zeros(3,params.Nsim);
            obj.QuatRefA_BodyArr = zeros(4,params.Nsim);

            % Reaction wheel
            obj.params.N_react = size(obj.params.gs_b_arr,2);
            fprintf('%d reaction wheels \n',obj.params.N_react);

            for r = 1:obj.params.N_react 
                % Set the irrelevant gg_b_arr and gt_b_arr
                Z = null(obj.params.gs_b_arr(:,r)');
                obj.params.gg_b_arr(:,r) = Z(:,1); % Arbitrary assignment of gt and gg
                obj.params.gt_b_arr(:,r) = Z(:,2); 
            end

            obj.params.Itot_Body = obj.params.J_SatBody_C; % TODO!! This is wrong
        end

        function [U,TorqueB_C_ResB] = GetControls(obj,PosRef,VelRef,QuatRefA_Body,OmegaRefA_Body)
            % Inputs: 
            % QuatRefA_Body = quaternion representing Fref wrt FA,
            % expressed in body frame
            % OmegaRefA_Body = reference angular vel (FR wrt FA), expressed
            % in body frame. Must match kinematically to QuatRefA_Body for
            % attitude control

            % Outputs:
            % TorqueB_C_ResB is for intermediate study, where we directly
            % control the torque on satellite body, relative to COM,
            % expressed in body frame. Set to empty if doing proper
            % actuator control

            % Here we compute standard attitude control terms, all in body
            % frame
            Xk = obj.statesArr(:,obj.SimCnt + 1);
            OmegaBA_Body = Xk(11:13);
            OmegeRefDot_Body = (OmegaRefA_Body - obj.OmegaRef_Prev)/obj.params.Ts;
            J_BC_OmegaBA_Body = obj.params.Itot_Body*OmegaBA_Body;
            Delta_W_Body = OmegaBA_Body - OmegaRefA_Body;
            if obj.SimCnt == 0 obj.Delta_W0_Body = Delta_W_Body; end;

            % quaternion error
            quat_BA_Body_k = Xk(7:10); % FB quat wrt FA
            O_BA_K = quat2dcm([quat_BA_Body_k(end);quat_BA_Body_k(1:3)]'); % [TODO: use my code and compare]
            O_RefA_K = quat2dcm([QuatRefA_Body(end);QuatRefA_Body(1:3)]');
            O_BR_K = O_BA_K*(O_RefA_K');

            if true % obtain thru DCM
                [err_quat_BR_Body_ScalerFirst] = dcm2quat(O_BR_K);
                err_quat_BR_Body = [err_quat_BR_Body_ScalerFirst(2:4)';err_quat_BR_Body_ScalerFirst(1)]; % convert to scaler last from Matlab
            else % direct quaternion [TODO: use direct quat and compare!]
            end
            err_quat_scaler_BR_Body = err_quat_BR_Body(end); err_quat_vec_BR_Body = err_quat_BR_Body(1:3);

            U = zeros(obj.params.N_react,1);
            TorqueB_C_ResB = [];

            switch obj.params.ControlScheme
                case 505 % atittude-rate control option 1 using P control 
                    % Implements the attitude-rate control option 1, essentially a
                    % nonlinear proportional controller of delta w = w_BA - w_ref
                    % resolved in satellite BODY frame
                    P = diag([1.0;1.0;1.0]);
                    TorqueB_C_ResB = cross(OmegaRefA_Body,J_BC_OmegaBA_Body) + obj.params.Itot_Body*OmegeRefDot_Body ...
                                    -P*Delta_W_Body;
                case 506 % 3D attitude control option 1 using PID
                    % compute the integral Z
                    Kp = 0.1; Ki = 0.00; Kd = 3.0;
                    obj.Body_int_506 = obj.Body_int_506 + SignJP(err_quat_scaler_BR_Body)*Kp*err_quat_vec_BR_Body*obj.params.Ts;
                    Z_body = obj.Body_int_506 + obj.params.Itot_Body*(Delta_W_Body - obj.Delta_W0_Body); 

                    TorqueB_C_ResB = cross(OmegaBA_Body,J_BC_OmegaBA_Body) + obj.params.Itot_Body*OmegeRefDot_Body ...
                                    -Kd*Delta_W_Body - Kp*SignJP(err_quat_scaler_BR_Body) *err_quat_vec_BR_Body - Kd*Ki*Z_body;
                otherwise
                    
            end
        end

        function Step(obj,PosRef,VelRef,QuatRefA_Body,OmegaRefA_Body)
            % 1 step in the current simulation

            % Compute controls
            [U,TorqueB_C_ResB] = obj.GetControls(PosRef,VelRef,QuatRefA_Body,OmegaRefA_Body);

            if obj.SimCnt == 0 % init UarrStore based on number of controls
                obj.UarrStore = zeros(length(U),obj.params.Nsim);
            end

            % propagate state forward 1 step
            odeFun = @(t,x) obj.satellite_dynamics(t, x,U, obj.planet, obj.params,TorqueB_C_ResB);
            tspan = [obj.SimT (obj.SimT+obj.params.Ts)]; 
            xCurrent = obj.statesArr(:,obj.SimCnt + 1);
            [T,X] = ode45(odeFun, tspan, xCurrent); 

           % Store memory and propagate counter
           obj.UarrStore(:,obj.SimCnt + 1) = U;
           obj.statesArr(:,obj.SimCnt + 2) = X(end,:)';
           obj.tArr(:,obj.SimCnt + 1) = obj.SimT;
           obj.OmegaRefA_BodyArr(:,obj.SimCnt + 1) = OmegaRefA_Body;
           obj.QuatRefA_BodyArr(:,obj.SimCnt + 1) = QuatRefA_Body;

           obj.SimCnt = obj.SimCnt + 1;
           obj.SimT = obj.SimCnt*obj.params.Ts;
           

           obj.OmegaRef_Prev = OmegaRefA_Body;
           if mod(obj.SimCnt,10) == 0 fprintf('%d\n',obj.SimCnt); end

        end        
       
    end % public

    methods (Static)
        function xdot = satellite_dynamics(t, x,U, planet, sat_params,TorqueB_C_ResB)
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
            external_torque = zeros(3,1); % Edit here for external torques [3;-2;1]*0.001;
            I_wdot = -cross(omega,sat_params.Itot_Body*omega) + external_torque;
            OmeReacDot = U(1:sat_params.N_react);

            % Get internal torque
            if ~isempty(TorqueB_C_ResB)
                I_wdot = I_wdot + TorqueB_C_ResB; 
            else                
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
            end

            omega_dot = sat_params.Itot_Body\I_wdot;
            xdot = [r_dot; v_dot; q_dot; omega_dot;OmeReacDot];
        end
    end % static
end

