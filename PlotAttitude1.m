figure(figN); figN = figN + 1; lm = []; kk = 1;
nRows = 3; nCols = 4;
RspCol = 'b'; CmdCol = 'r--'; mLineWidth = 1.5;

quats = statesArr(7:10,:)';
euls = quat2eul([quats(:, 4), quats(:, 1:3)], 'ZYX'); % convert to scalar first then to euler angles

pRadsArr= statesArr(11,:); qRadsArr= statesArr(12,:); rRadsArr= statesArr(13,:);
w_BA_B_RadsMagArr = sqrt(pRadsArr.^2 + qRadsArr.^2 + rRadsArr.^2);
w_CmdA_B_RadsMagArr = sqrt(OmegaCmdA_BodyArr(1,:).^2 + OmegaCmdA_BodyArr(2,:).^2 + OmegaCmdA_BodyArr(3,:).^2);

% Angular velocity
lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,pRadsArr*180/pi,RspCol,tArr,OmegaCmdA_BodyArr(1,:)*180/pi,CmdCol, 'LineWidth', mLineWidth);
ylabel('p (deg/s)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,qRadsArr*180/pi,RspCol,tArr,OmegaCmdA_BodyArr(2,:)*180/pi,CmdCol, 'LineWidth', mLineWidth);
ylabel('q (deg/s)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,rRadsArr*180/pi,RspCol,tArr,OmegaCmdA_BodyArr(3,:)*180/pi,CmdCol, 'LineWidth', mLineWidth);
ylabel('r (deg/s)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,w_BA_B_RadsMagArr*180/pi,RspCol, 'LineWidth', mLineWidth);
plot(tArr,w_CmdA_B_RadsMagArr*180/pi,CmdCol, 'LineWidth', mLineWidth);
ylabel('$w_{BA}^B$ mag (deg/s)','Interpreter','latex');


lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,1),RspCol,tArr,QuatRefA_A_Arr(1,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q1');


lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,2),RspCol,tArr,QuatRefA_A_Arr(2,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q2');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,3),RspCol,tArr,QuatRefA_A_Arr(3,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q3');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,4),RspCol,tArr,QuatRefA_A_Arr(4,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q4');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,Err_BR_AngRadArr*180/pi,RspCol, 'LineWidth', mLineWidth);
ylabel('Angle error [deg]');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,euls*180/pi, 'LineWidth', mLineWidth)
legend("$\psi$", "$\theta$", "$\phi$", 'Interpreter', 'latex');
title('Angular Position');
ylabel('(deg)');

% Plot reaction wheels
lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
leg = [];
title('Reaction wheel angular velocity');
for r = 1:SatellitePlot.params.N_react
    leg(r) = plot(tArr,statesArr(13+r,:)/(2*pi)*60, 'LineWidth', mLineWidth,'DisplayName',num2str(r));   
end
% yline(6500,'r--'); yline(-6500,'r--'); % NanoAvio Max RPM
ylabel('(RPM)');
xlabel('Time (s)');
legend(leg, 'Interpreter', 'latex');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
leg = [];
title('Reaction wheel angular acceleration');
for r = 1:SatellitePlot.params.N_react
    leg(r) = plot(tArr,UarrStore(r,:)/(2*pi)*60^2, 'LineWidth', mLineWidth,'DisplayName',num2str(r));   
end
MaxTorque = 3.2/1000; % NanoAvio Max torque
MaxReacWheelAngAccRads2 = MaxTorque/sat_params.I_ws;
% yline(MaxReacWheelAngAccRads2/(2*pi)*60^2,'r--'); yline(-MaxReacWheelAngAccRads2/(2*pi)*60^2,'r--'); 
ylabel('(RPM2)');
xlabel('Time (s)');
legend(leg, 'Interpreter', 'latex');

linkaxes(lm,'x');