% ALL checked!
%% define desired refA 
% Delta rotation vector and angles wrt latest desired reference orientation 
DeltaRotVec_CurRef = [1 1 1;
               1 1 -1;
               1 -1 1;
               1 -1 -1;
               -1 1 1;
               -1 1 -1;
               -1 -1 1;
               -1 -1 -1;];
DeltaRotAngDeg = [45 90 135 180 180 180 180 -120]';

if false
    DeltaRotVec_CurRef = [0 0 1;0 1 0;1 0 0];
    DeltaRotAngDeg = [90;-30;45];
end

nChg = length(DeltaRotAngDeg);

if k == 0
    QuatRefA_A_Desired_Scal = Q0_A(4);
    QuatRefA_A_Desired_Vec = Q0_A(1:3);
    ChgIdx = 1;
else
    O_RefA_K = quat2dcm([QuatRefA_A_Desired_Scal;QuatRefA_A_Desired_Vec]'); % [TODO: use my code and compare]

    % Define change in desired quat ref 
    ChangePeriodCnt = 100;
    if mod(k-1,ChangePeriodCnt) == 0
        if ChgIdx <= nChg
            fprintf('Chg %d\n',k);
            [DeltaQuat_Scaler,DeltaQuat_Vec_A] = RotVectAngle2Quat(O_RefA_K'*DeltaRotVec_CurRef(ChgIdx,:)',DeltaRotAngDeg(ChgIdx)*pi/180);
            [QuatRefA_A_Desired_Scal,QuatRefA_A_Desired_Vec] = QuatHamDot(DeltaQuat_Scaler,DeltaQuat_Vec_A,...
                QuatRefA_A_Desired_Scal,QuatRefA_A_Desired_Vec); % q_DesRefNewA = q_DesRefNew_DesRefOld (dot) q_DesRefOld
            ChgIdx = ChgIdx + 1;
        end
    end
end

%% Slew QuatRefA to QuatRefA_Body_Desired
SlewRateLimitDegs = 1.0;

[qiR_scal,qiR_vec,RotVec_21,AngRad_i1] = InterpolateQuat_PhysicalAngularRateLimit(...
    QuatRefA_A(4),QuatRefA_A(1:3),QuatRefA_A_Desired_Scal,QuatRefA_A_Desired_Vec,SlewRateLimitDegs*pi/180*Ts);

QuatRefA_A = [qiR_vec;qiR_scal];