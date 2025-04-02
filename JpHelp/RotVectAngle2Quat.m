function [Quat_Scaler,Quat_Vec] = RotVectAngle2Quat(RotVec,AngRad)
    % Convert rotation's vector RotVec and angle AngRad to quaternion's
    % Quat_Scaler and Quat_Vec
    Quat_Scaler = cos(AngRad/2);
    Quat_Vec = sin(AngRad/2)*RotVec;
end

