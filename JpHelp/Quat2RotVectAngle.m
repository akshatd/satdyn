function [RotVec,AngRad] = Quat2RotVectAngle(Quat_Scaler,Quat_Vec)
    % Convert quat's Quat_Scaler and Quat_Vec to rotation's vector RotVec and angle AngRad
    cos_ = Quat_Scaler;
    sin_ = norm(Quat_Vec);

    AngRad = atan2(sin_,cos_)*2;
    if abs(sin_) > 0.00001
        RotVec = Quat_Vec/sin_;
    else
        RotVec = [1;0;0]; % arbitrary
    end

end

