function [YawR,TetR,RollR] = DCM_to_Euler321(O_BA)
    % convert DCM to 321 yaw pitch roll euler sequence 
    % DCM = O_BA, the orientation matrix to change perspective X_B =
    % O_BA*X_A, where X_B is vector X resolved in frame B

    YawR = atan2(O_BA(1,2),O_BA(1,1));
    TetR =-asin(O_BA(1,3));
    RollR = atan2(O_BA(2,3),O_BA(3,3));
end

