% Define quaternion
nta_BA = [1;0;0.0]; tetd_BA = 90;
nta_BA = [0.2;0.45;-0.9]; tetd_BA = 40;

nta_BA = nta_BA/norm(nta_BA);

scal_BA = cosd(tetd_BA/2);
vect_BA = sind(tetd_BA/2)*nta_BA;

%% quat to DCM
DCM_mine = quat2dcm_jp(scal_BA,vect_BA);
DCM_Matlab = quat2dcm([scal_BA;vect_BA]' ); % same as MATLAB

%% DCM to quad
[scal_BA_mine,vec_BA_mine] = dcm2quat_jp(DCM_mine);
[quat_BA_Matlab_ScalerFirst] = dcm2quat(DCM_Matlab); % same as MATLAB

%% DCM to euler
[YawR_mine,TetR_mine,RollR_mine] = DCM_to_Euler321(DCM_mine);
[YawR_Matlab TetR_Matlab RollR_Matlab] = dcm2angle(DCM_Matlab,'ZYX'); % same as MATLAB

%% Euler to DCM
DCM_mine2 = Euler321_to_dcm(YawR_mine,TetR_mine,RollR_mine);
DCM_Matlab2 = angle2dcm(YawR_Matlab ,TetR_Matlab ,RollR_Matlab,'ZYX');  % same as MATLAB

%% Euler to quat
[scal_BA_mine2,vec_BA_mine2] = dcm2quat_jp(Euler321_to_dcm(YawR_mine,TetR_mine,RollR_mine));
quat_BA_Matlab_ScalerFirst2 = eul2quat([YawR_Matlab ,TetR_Matlab ,RollR_Matlab], 'ZYX');  % same as MATLAB

%% quat to Euler
[YawR_mine,TetR_mine,RollR_mine] = DCM_to_Euler321(quat2dcm_jp(scal_BA_mine2,vec_BA_mine2));
euls_Matlab = quat2eul(quat_BA_Matlab_ScalerFirst2, 'ZYX'); % euls_Matlab = YPR , % same as MATLAB