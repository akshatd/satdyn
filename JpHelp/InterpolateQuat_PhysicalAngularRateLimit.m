function [qiR_scal,qiR_vec,RotVec_21,AngRad_i1] = InterpolateQuat_PhysicalAngularRateLimit(q1R_Scal,q1R_vec,q2R_Scal,q2R_vec,MaxPhysicalAngChangeLimit_Rad)
    % interpolate from q1R to q2R (defines rotation wrt common frame FR)
    % with a maximum physical angle change limit MaxPhysicalAngChangeLimit_Rad
    % returns interpolated quat (qiR_scal,qiR_vec) and rotation axis
    % RotVec_21 and actual rotated angle AngRad_i1

    % Refer to my 'Quaternions Detailed: Rot quat interp with physical
    % angle given'

    [q21_scal,q21_vec] = QuatHamDot(q2R_Scal,q2R_vec,q1R_Scal,-q1R_vec); % q2R dot q1R inv = q21

    % Get the rotational axis and angle to go from q1R to q2R
    [RotVec_21,AngRad_21] = Quat2RotVectAngle(q21_scal,q21_vec);

    % limit the actual rotation angle
    AngRad_i1 = median([-MaxPhysicalAngChangeLimit_Rad AngRad_21 MaxPhysicalAngChangeLimit_Rad]);
    
    % Get qi1
    [qi1_Scaler,qi1_Vec] = RotVectAngle2Quat(RotVec_21,AngRad_i1);

    % Get interpolated qiR
    [qiR_scal,qiR_vec] = QuatHamDot(qi1_Scaler,qi1_Vec,q1R_Scal,q1R_vec);
end

