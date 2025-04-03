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

        TorqueB_C_ResB_Prev
        U_Prev

        OmegaCmdA_BodyArr   % store OmegaCmdA_Body
        QuatRefA_A_Arr   % store QuatRefA_A

        OmegaRefA_Body_Prev = zeros(3,1);

        Body_int_506 = zeros(3,1); % integrator vector state for option 506
        Delta_W0_Body % store initial delta W

        UarrGuess_Prev_507  % previous iteration's optimal guess

        W_BA_Cmd_Body_Prev_508  % previous W_BA_Cmd_Body

        UarrGuess_Prev_509  % previous iteration's optimal guess

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
            % gs_b_arr_MinNorm = [gs_b_arr]'*inv(gs_b_arr*gs_b_arr')


            % init states storage for entire simulation
          

            obj.params = params;
            obj.planet = planet;

            % Define initial state vector
            states0 = [R0;V0;Q0;W0;W_React0arr];

            obj.statesArr = zeros(length(states0),params.Nsim + 1);
            obj.statesArr(:,1) = states0;
            obj.tArr = zeros(1,params.Nsim);
            obj.OmegaCmdA_BodyArr = zeros(3,params.Nsim);
            obj.QuatRefA_A_Arr = zeros(4,params.Nsim);

            % Reaction wheel
            obj.params.N_react = size(obj.params.gs_b_arr,2);
            fprintf('%d reaction wheels \n',obj.params.N_react);

            % for r = 1:obj.params.N_react 
            %     % Set the irrelevant gg_b_arr and gt_b_arr
            %     Z = null(obj.params.gs_b_arr(:,r)');
            %     obj.params.gg_b_arr(:,r) = Z(:,1); % Arbitrary assignment of gt and gg
            %     obj.params.gt_b_arr(:,r) = Z(:,2); 
            % end
            obj.params.gs_b_arr_MinNorm = [obj.params.gs_b_arr]'*inv(obj.params.gs_b_arr*obj.params.gs_b_arr');

            obj.params.Itot_Body = obj.params.J_SatBody_C; % NOTE!! Assumes J_SatBody_C includes ALL N_react reaction wheel's contribution too for simplicity
            obj.params.Itot_inv_Body = inv(obj.params.Itot_Body);
        end

        function [U,TorqueB_C_ResB] = GetControls(obj,PosRef,VelRef,QuatRefA_A,OmegaRefA_A)
            % Inputs: 
            % QuatRefA_A = quaternion representing Fref wrt FA,
            % expressed in inertial frame
            % OmegaRefA_A = reference angular vel (FR wrt FA), expressed
            % in inertial frame. Must match kinematically to QuatRefA_A for
            % attitude control

            % Outputs:
            % TorqueB_C_ResB is for intermediate study, where we directly
            % control the torque on satellite body, relative to COM,
            % expressed in body frame. Set to empty if doing proper
            % actuator control

            W_mag_rads_limit = 1*pi/180; Umag_rads2_limit = 88.3087;% NanoAvio 1*pi/180;

            if obj.SimCnt == 0
                obj.TorqueB_C_ResB_Prev = zeros(3,1);
                obj.U_Prev = zeros(obj.params.N_react,1);
            end

            Xk = obj.statesArr(:,obj.SimCnt + 1);
            quat_BA_Body_k = Xk(7:10); % FB quat wrt FA
            O_BA_K = quat2dcm([quat_BA_Body_k(end);quat_BA_Body_k(1:3)]'); % [TODO: use my code and compare]

            OmegaRefA_Body = O_BA_K*OmegaRefA_A; % TODO: use quat to change perspective
            if obj.SimCnt == 0 obj.OmegaRefA_Body_Prev = OmegaRefA_Body; end;

            % Here we compute standard attitude control terms, all in body
            % frame            
            OmegaBA_Body = Xk(11:13);
            OmegeRefA_Dot_Body = (OmegaRefA_Body - obj.OmegaRefA_Body_Prev)/obj.params.Ts;
            J_BC_OmegaBA_Body = obj.params.Itot_Body*OmegaBA_Body;
            Delta_W_Body = OmegaBA_Body - O_BA_K*OmegaRefA_A; % TODO: use quat to change perspective
            if obj.SimCnt == 0 obj.Delta_W0_Body = Delta_W_Body; end;

            % quaternion error                    
            O_RefA_K = quat2dcm([QuatRefA_A(end);QuatRefA_A(1:3)]');
            O_BR_K = O_BA_K*(O_RefA_K');

            if true % obtain thru DCM
                [err_quat_BR_Body_ScalerFirst] = dcm2quat(O_BR_K);
                err_quat_BR_Body = [err_quat_BR_Body_ScalerFirst(2:4)';err_quat_BR_Body_ScalerFirst(1)]; % convert to scaler last from Matlab
            else % direct quaternion [TODO: use direct quat and compare!]
            end
            err_quat_scaler_BR_Body = err_quat_BR_Body(end); err_quat_vec_BR_Body = err_quat_BR_Body(1:3);

            % Predict ahead the ref quat using OmegaRef
            ell = 10; % MPC horizon
            xRefArr = ODE_RK4(@(X,U,params) SatelliteClass.fRef_cont(X),[QuatRefA_A;O_RefA_K*OmegaRefA_A],ell,zeros(1,ell),obj.params.Ts,[]);
            xRefArr = xRefArr(:,2:end); % Xref = [QuatRefA_A;OmegaRefA_Ref] propagated with constant OmegaRefA_Ref

            % Reaction wheel matrices
            omega_reac = Xk(14:(13+obj.params.N_react));
            hW = obj.params.I_ws*omega_reac;
            TorqueReacFF_ResB = cross(OmegaBA_Body,obj.params.gs_b_arr*hW);


            U = zeros(obj.params.N_react,1);
            TorqueB_C_ResB = [];

            % General control tuning
            Tau_timeConstant_Secs = 30; Zeta = 0.9;

            switch obj.params.ControlScheme
                case 505 % atittude-rate control option 1 using P control                    
                    TorqueB_C_ResB = SatelliteClass.GetAttitudeRateControl_Option1_Pctrl(OmegaRefA_Body,...
                        OmegeRefA_Dot_Body,J_BC_OmegaBA_Body,OmegaBA_Body,obj.params);

                    Tau = TorqueB_C_ResB + TorqueReacFF_ResB;
                    U = -obj.params.gs_b_arr_MinNorm *Tau/obj.params.I_ws;
                    U(U > Umag_rads2_limit) = Umag_rads2_limit;
                    U(U < -Umag_rads2_limit) = -Umag_rads2_limit; %if norm(U) > Umag_rads2_limit U = U*Umag_rads2_limit/norm(U);end
                    
                    TorqueB_C_ResB = [];
                    
                    OmegaCmdA_Body = OmegaRefA_Body;
                case 506 % 3D attitude control option 1 using PID
                    % compute the integral Z
                    Kd = 2*obj.params.J_SatBody_C_Mean/Tau_timeConstant_Secs;
                    Kp = (Kd/Zeta)^2/2/obj.params.J_SatBody_C_Mean;
                    Ki = 0.00; %Kp = 0.05; Kd = 1.0;
                    obj.Body_int_506 = obj.Body_int_506 + SignJP(err_quat_scaler_BR_Body)*Kp*err_quat_vec_BR_Body*obj.params.Ts;
                    Z_body = obj.Body_int_506 + obj.params.Itot_Body*(Delta_W_Body - obj.Delta_W0_Body); 

                    TorqueB_C_ResB = cross(OmegaBA_Body,J_BC_OmegaBA_Body) + obj.params.Itot_Body*OmegeRefA_Dot_Body ...
                                    -Kd*Delta_W_Body - Kp*SignJP(err_quat_scaler_BR_Body) *err_quat_vec_BR_Body - Kd*Ki*Z_body;

                    Tau = TorqueB_C_ResB + TorqueReacFF_ResB;
                    U = -obj.params.gs_b_arr_MinNorm *Tau/obj.params.I_ws;
                    U(U > Umag_rads2_limit) = Umag_rads2_limit;
                    U(U < -Umag_rads2_limit) = -Umag_rads2_limit; %if norm(U) > Umag_rads2_limit U = U*Umag_rads2_limit/norm(U);end
                    
                    TorqueB_C_ResB = [];

                    OmegaCmdA_Body = OmegaRefA_Body;

                case 507 % 3D attitude control option 2 using basic nonlinear MPC
                    options = optimoptions('fmincon','Display','iter','MaxIterations',50,'Algorithm','sqp');
                    Xk_MPC2 = [quat_BA_Body_k;OmegaBA_Body];
                    
                    Q_err_quat_vec_BR_Body = 0.1*eye(3);
                    Q_Delta_W_Body = 0.5*eye(3);
                    R_Delta_U  = 0.1*eye(3);

                    obj_ = @(Uarr) SatelliteClass.obj_fun_MPC_Option2(obj.params,Uarr,Xk_MPC2,obj.TorqueB_C_ResB_Prev,ell,xRefArr,Q_err_quat_vec_BR_Body,Q_Delta_W_Body,R_Delta_U);
                    NONLCON = @(Uarr) SatelliteClass.constr_fun_MPC_Option2(Uarr,Xk_MPC2,ell,W_mag_rads_limit,obj.params);
                    % fmincon(J,W0,A1,B1,A2,B2,LB_W,UB_W,NONLCON), min(W) J(W) s.t. A1*W  <= B1, A2*W  = B2
                    % LB_W <= W <= UB_W. [g,h] = NONLCON(W) where g(W) <= 0, h(W) = 0 and 
                    if obj.SimCnt ==0 obj.UarrGuess_Prev_507 = ones(3*ell,1)*0.01;end
                    UB = ones(3*ell,1)*inf;
                    LB = -UB;
                    Uarr_opt = fmincon(obj_,[obj.UarrGuess_Prev_507(4:end); obj.UarrGuess_Prev_507(end-2:end)],[],[],[],[],LB,UB,NONLCON,options);
                    obj.UarrGuess_Prev_507 = Uarr_opt;
                    
                    TorqueB_C_ResB = Uarr_opt(1:3);
                    OmegaCmdA_Body = OmegaRefA_Body;

                case 508 % 3D attitude control option 3 (P-outer, W-inner)
                    error('Not complete\n')
                    Kp = 0.1; 

                    % Compute W_BA_Cmd_Body
                    delta_WrefA_Body = -Kp*SignJP(err_quat_scaler_BR_Body)*err_quat_vec_BR_Body;
                    W_BA_Cmd_Body = delta_WrefA_Body + OmegaRefA_Body;
                    Wcmd_mag_rads = norm(W_BA_Cmd_Body);
                    if Wcmd_mag_rads > W_mag_rads_limit
                        W_BA_Cmd_Body = W_BA_Cmd_Body*W_mag_rads_limit/Wcmd_mag_rads;
                    end
                    if obj.SimCnt == 0 obj.W_BA_Cmd_Body_Prev_508 = W_BA_Cmd_Body;end

                    OmegeCmdA_Dot_Body = (W_BA_Cmd_Body - obj.W_BA_Cmd_Body_Prev_508)/obj.params.Ts;

                    % inner-loop
                    TorqueB_C_ResB = SatelliteClass.GetAttitudeRateControl_Option1_Pctrl(W_BA_Cmd_Body,...
                        OmegeCmdA_Dot_Body,J_BC_OmegaBA_Body,OmegaBA_Body,obj.params);

                    obj.W_BA_Cmd_Body_Prev_508 = W_BA_Cmd_Body;
                case 509 % 3D attitude control option 4 using reaction wheels nonlinear MPC
                    options = optimoptions('fmincon','Display','iter','MaxIterations',50,'Algorithm','active-set');
                    Xk_MPC4 = [quat_BA_Body_k;OmegaBA_Body;omega_reac];
                    
                    Q_err_quat_vec_BR_Body = 1.5*eye(3);
                    Q_Delta_W_Body = 300.0*eye(3);
                    R_Delta_U  = 0.01*eye(obj.params.N_react);
                    R_U = 0.005*eye(obj.params.N_react);

                    if obj.SimCnt ==0 obj.UarrGuess_Prev_509 = ones(obj.params.N_react*ell,1)*0.01;end

                    Ukm1_MPC4 = obj.UarrGuess_Prev_509(1:obj.params.N_react);

                    obj_ = @(Uarr) SatelliteClass.obj_fun_MPC_Option4(obj.params,Uarr,Xk_MPC4,Ukm1_MPC4,ell,xRefArr,Q_err_quat_vec_BR_Body,Q_Delta_W_Body,R_Delta_U,R_U);
                    NONLCON = @(Uarr) SatelliteClass.constr_fun_MPC_Option4(Uarr,Xk_MPC4,ell,W_mag_rads_limit,obj.params);
                    % fmincon(J,W0,A1,B1,A2,B2,LB_W,UB_W,NONLCON), min(W) J(W) s.t. A1*W  <= B1, A2*W  = B2
                    % LB_W <= W <= UB_W. [g,h] = NONLCON(W) where g(W) <= 0, h(W) = 0 and 
                    
                    UB = ones(obj.params.N_react*ell,1)*Umag_rads2_limit;
                    LB = -UB;
                    Uarr_opt = fmincon(obj_,[obj.UarrGuess_Prev_509((obj.params.N_react + 1):end); obj.UarrGuess_Prev_509((end-obj.params.N_react + 1):end)],[],[],[],[],LB,UB,NONLCON,options);
                    obj.UarrGuess_Prev_509 = Uarr_opt;
                    
                    U = Uarr_opt(1:obj.params.N_react); TorqueB_C_ResB = [];

                    OmegaCmdA_Body = OmegaRefA_Body;

                otherwise
            end

            obj.OmegaRefA_Body_Prev = OmegaRefA_Body;
            obj.OmegaCmdA_BodyArr(:,obj.SimCnt + 1) = OmegaCmdA_Body;
        end

        function Step(obj,PosRef,VelRef,QuatRefA_A,OmegaRefA_A)
            % 1 step in the current simulation

            % Compute controls
            [U,TorqueB_C_ResB] = obj.GetControls(PosRef,VelRef,QuatRefA_A,OmegaRefA_A);

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
           
           obj.QuatRefA_A_Arr(:,obj.SimCnt + 1) = QuatRefA_A;
           obj.TorqueB_C_ResB_Prev = TorqueB_C_ResB;
           obj.U_Prev = U;

           obj.SimCnt = obj.SimCnt + 1;
           obj.SimT = obj.SimCnt*obj.params.Ts;
           

           
           if mod(obj.SimCnt,10) == 0 fprintf('%d\n',obj.SimCnt); end

        end        
    end % public

    methods (Static)
        
        function Xd = fRef_cont(Xref)
            % Xref = [QuatRefA_A;OmegaRefA_Ref]
            % just calculate quat dot, assuming OmegaRefA_Ref dot = 0
            QuatRefA_A = Xref(1:4);
            OmegaRefA_Ref = Xref(5:7);

            q_dot = 0.5 * [
	            QuatRefA_A(4), -QuatRefA_A(3), QuatRefA_A(2);
	            QuatRefA_A(3), QuatRefA_A(4), -QuatRefA_A(1);
	            -QuatRefA_A(2), QuatRefA_A(1), QuatRefA_A(4);
	            -QuatRefA_A(1), -QuatRefA_A(2), -QuatRefA_A(3)
	            ] * OmegaRefA_Ref;
            Xd = [q_dot;zeros(3,1)];
        end

        % Start ---------- 3D attitude control option 4 (reaction wheel Nonlin MPC -----------------

        % Define continous-time dynamics
        function Xd = f_cont_MPC_Option4(X_MPC4,U_MPC4,params)
            q_BA_B = X_MPC4(1:4); % q_BA_(B or A)
            W_BA_B = X_MPC4(5:7);     % W_BA_B
            omega_reac = X_MPC4(8:(7 + params.N_react)); % reaction wheels

            q_dot = 0.5 * [
	            q_BA_B(4), -q_BA_B(3), q_BA_B(2);
	            q_BA_B(3), q_BA_B(4), -q_BA_B(1);
	            -q_BA_B(2), q_BA_B(1), q_BA_B(4);
	            -q_BA_B(1), -q_BA_B(2), -q_BA_B(3)
	            ] * W_BA_B;

            hW = params.I_ws*omega_reac;
            TorqueReacFF_ResB = cross(W_BA_B,params.gs_b_arr*hW);

            W_BA_B_dot = params.Itot_inv_Body*(-cross(W_BA_B,params.Itot_Body*W_BA_B) - TorqueReacFF_ResB...
                                                -params.gs_b_arr*params.I_ws*U_MPC4);
            Xd = [q_dot; W_BA_B_dot;U_MPC4];            
        end
        % Define objective function
        function J = obj_fun_MPC_Option4(params,Uarr,Xk_MPC4,Ukm1_MPC4,ell,Xref_Arr,Q_err_quat_vec_BR_Body,Q_Delta_W_Body,R_Delta_U,R_U)
            Uarr = reshape(Uarr,params.N_react,ell);

            QuatRefA_A_Arr = Xref_Arr(1:4,:);
            OmegaRefA_Ref_Arr = Xref_Arr(5:7,:);

            X_Arr = ODE_RK4(@SatelliteClass.f_cont_MPC_Option4,Xk_MPC4,ell,Uarr,params.Ts,params);
            X_Arr = X_Arr(:,2:end);           

            % delta w error in body frame (assume to be approx ref frame at every prediction step)
            OmegaBA_Body_Arr = X_Arr(5:7,:);
            Delta_W_Body_kArr = OmegaBA_Body_Arr - OmegaRefA_Ref_Arr;

            % quaternion error
            quat_BA_Body_kArr = X_Arr(1:4,:); % FB quat wrt FA
            err_quat_BR_Body_kArr = zeros(4,ell);
            for i = 1:ell
                QuatRefA_A_k = QuatRefA_A_Arr(:,i);
                quat_BA_Body_k = quat_BA_Body_kArr(:,i);
                O_BA_K = quat2dcm([quat_BA_Body_k(end);quat_BA_Body_k(1:3)]'); % [TODO: use my code and compare]
                O_RefA_K = quat2dcm([QuatRefA_A_k(end);QuatRefA_A_k(1:3)]');
                O_BR_K = O_BA_K*(O_RefA_K');
    
                if true % obtain thru DCM
                    [err_quat_BR_Body_ScalerFirst] = dcm2quat(O_BR_K);
                    err_quat_BR_Body = [err_quat_BR_Body_ScalerFirst(2:4)';err_quat_BR_Body_ScalerFirst(1)]; % convert to scaler last from Matlab
                else % direct quaternion [TODO: use direct quat and compare!]
                end
                err_quat_BR_Body_kArr(:,i) = err_quat_BR_Body;

            end

            % Control rate
            Delta_U_kArr = diff([Ukm1_MPC4 Uarr],1,2);

            % Compute cost function
            J = 0;
            for i = 0:(ell-1)
                J = J + err_quat_BR_Body_kArr(1:3,i+1)'*Q_err_quat_vec_BR_Body*err_quat_BR_Body_kArr(1:3,i+1) + ...
                        Delta_W_Body_kArr(:,i+1)'*Q_Delta_W_Body*Delta_W_Body_kArr(:,i+1) + ...
                        Delta_U_kArr(:,i+1)'*R_Delta_U*Delta_U_kArr(:,i+1) + ...
                        Uarr(:,i+1)'*R_U*Uarr(:,i+1);
            end            
        end

        % Define constraint function
        function [g,h] = constr_fun_MPC_Option4(Uarr,Xk_MPC4,ell,W_mag_rads_limit,params)
            % Note: g(Uarr) <=0, h(Uarr) = 0
            Uarr = reshape(Uarr,params.N_react,ell);

            X_Arr = ODE_RK4(@SatelliteClass.f_cont_MPC_Option4,Xk_MPC4,ell,Uarr,params.Ts,params);
            X_Arr = X_Arr(:,2:end);    
            

            pRadsArr= X_Arr(5,:); qRadsArr=  X_Arr(6,:); rRadsArr=  X_Arr(7,:);
            w_BA_B_RadsMagArr = sqrt(pRadsArr.^2 + qRadsArr.^2 + rRadsArr.^2);

            g = w_BA_B_RadsMagArr -W_mag_rads_limit;            
            h = [];
        end
        % End ---------- 3D attitude control option 4 (reaction wheel Nonlin MPC -----------------

        % functions for 3D attitude control option 2 (basic nonlinear MPC)
        % Define continous-time dynamics
        function Xd = f_cont_MPC_Option2(X_MPC2,U_MPC2,params)
            q_BA_B = X_MPC2(1:4); % q_BA_(B or A)
            W_BA_B = X_MPC2(5:7);     % W_BA_B
            q_dot = 0.5 * [
	            q_BA_B(4), -q_BA_B(3), q_BA_B(2);
	            q_BA_B(3), q_BA_B(4), -q_BA_B(1);
	            -q_BA_B(2), q_BA_B(1), q_BA_B(4);
	            -q_BA_B(1), -q_BA_B(2), -q_BA_B(3)
	            ] * W_BA_B;

            Xd = [q_dot;
            params.Itot_inv_Body*(U_MPC2 - cross(W_BA_B,params.Itot_Body*W_BA_B))];            
        end
        % Define objective function
        function J = obj_fun_MPC_Option2(params,Uarr,Xk_MPC2,Ukm1_MPC2,ell,Xref_Arr,Q_err_quat_vec_BR_Body,Q_Delta_W_Body,R_Delta_U)
            Uarr = reshape(Uarr,3,ell);

            QuatRefA_A_Arr = Xref_Arr(1:4,:);
            OmegaRefA_Ref_Arr = Xref_Arr(5:7,:);

            X_Arr = ODE_RK4(@SatelliteClass.f_cont_MPC_Option2,Xk_MPC2,ell,Uarr,params.Ts,params);
            X_Arr = X_Arr(:,2:end);           

            % delta w error in body frame (assume to be approx ref frame at every prediction step)
            OmegaBA_Body_Arr = X_Arr(5:7,:);
            Delta_W_Body_kArr = OmegaBA_Body_Arr - OmegaRefA_Ref_Arr;

            % quaternion error
            quat_BA_Body_kArr = X_Arr(1:4,:); % FB quat wrt FA
            err_quat_BR_Body_kArr = zeros(4,ell);
            for i = 1:ell
                QuatRefA_A_k = QuatRefA_A_Arr(:,i);
                quat_BA_Body_k = quat_BA_Body_kArr(:,i);
                O_BA_K = quat2dcm([quat_BA_Body_k(end);quat_BA_Body_k(1:3)]'); % [TODO: use my code and compare]
                O_RefA_K = quat2dcm([QuatRefA_A_k(end);QuatRefA_A_k(1:3)]');
                O_BR_K = O_BA_K*(O_RefA_K');
    
                if true % obtain thru DCM
                    [err_quat_BR_Body_ScalerFirst] = dcm2quat(O_BR_K);
                    err_quat_BR_Body = [err_quat_BR_Body_ScalerFirst(2:4)';err_quat_BR_Body_ScalerFirst(1)]; % convert to scaler last from Matlab
                else % direct quaternion [TODO: use direct quat and compare!]
                end
                err_quat_BR_Body_kArr(:,i) = err_quat_BR_Body;

            end

            % Control rate
            Delta_U_kArr = diff([Ukm1_MPC2 Uarr],1,2);

            % Compute cost function
            J = 0;
            for i = 0:(ell-1)
                J = J + err_quat_BR_Body_kArr(1:3,i+1)'*Q_err_quat_vec_BR_Body*err_quat_BR_Body_kArr(1:3,i+1) + ...
                        Delta_W_Body_kArr(:,i+1)'*Q_Delta_W_Body*Delta_W_Body_kArr(:,i+1) + ...
                        Delta_U_kArr(:,i+1)'*R_Delta_U*Delta_U_kArr(:,i+1);
            end            
        end

        % Define constraint function
        function [g,h] = constr_fun_MPC_Option2(Uarr,Xk_MPC2,ell,W_mag_rads_limit,params)
            % Note: g(Uarr) <=0, h(Uarr) = 0
            Uarr = reshape(Uarr,3,ell);

            X_Arr = ODE_RK4(@SatelliteClass.f_cont_MPC_Option2,Xk_MPC2,ell,Uarr,params.Ts,params);
            X_Arr = X_Arr(:,2:end);     

            pRadsArr= X_Arr(5,:); qRadsArr=  X_Arr(6,:); rRadsArr=  X_Arr(7,:);
            w_BA_B_RadsMagArr = sqrt(pRadsArr.^2 + qRadsArr.^2 + rRadsArr.^2);

            g = w_BA_B_RadsMagArr -W_mag_rads_limit;            
            h = [];
        end

        function TorqueB_C_ResB = GetAttitudeRateControl_Option1_Pctrl(OmegaCmdA_Body,OmegeCmdA_Dot_Body,J_BC_OmegaBA_Body,OmegaBA_Body,params)
            % atittude-rate control option 1 using P control 
            % Implements the attitude-rate control option 1, essentially a
            % nonlinear proportional controller of delta w = w_BA - w_CmdA
            % resolved in satellite BODY frame
            % Note that the command is OmegaCmdA
            Delta_W_Body = OmegaBA_Body - OmegaCmdA_Body;
            Pgain_Matrix = diag([1.0;1.0;1.0])*0.000001;
            TorqueB_C_ResB = cross(OmegaCmdA_Body,J_BC_OmegaBA_Body) + params.Itot_Body*OmegeCmdA_Dot_Body ...
                            -Pgain_Matrix*Delta_W_Body;
        end

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
            external_torque = zeros(3,1); % Edit here for external torques 
            % external_torque = [3;-2;1]*0.001;
            I_wdot = -cross(omega,sat_params.Itot_Body*omega) + external_torque;
            OmeReacDot = U(1:sat_params.N_react);

            % Get internal torque
            if ~isempty(TorqueB_C_ResB)
                I_wdot = I_wdot + TorqueB_C_ResB; 
            else     
                % vectorized reaction wheel's torque
                hW = sat_params.I_ws*omega_reac;
                I_wdot = I_wdot - cross(omega,sat_params.gs_b_arr*hW) - sat_params.gs_b_arr*sat_params.I_ws*OmeReacDot;

                % for r = 1:sat_params.N_react % Sum impact of each reaction wheel
                %     % Compute omega_G = [ws;wt;wg] = angular velocity of
                %     % satellite relative to inertial, resolved in gimbal frame
                % 
                %     % get O_BG
                %     O_BG = [sat_params.gs_b_arr(:,r) sat_params.gt_b_arr(:,r) sat_params.gg_b_arr(:,r)];
                %     omega_G = O_BG'*omega; 
                % 
                %     I_wdot = I_wdot - sat_params.I_ws*OmeReacDot(r)*sat_params.gs_b_arr(:,r) ...
                %                     - sat_params.I_ws*omega_reac(r)*omega_G(3)*sat_params.gt_b_arr(:,r) ...
                %                     + sat_params.I_ws*omega_reac(r)*omega_G(2)*sat_params.gg_b_arr(:,r);
                % 
                % end
            end

            omega_dot = sat_params.Itot_inv_Body*I_wdot;
            xdot = [r_dot; v_dot; q_dot; omega_dot;OmeReacDot];
        end
    end % static
end

