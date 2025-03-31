figure(figN); figN = figN + 1; lm = []; kk = 1;
nRows = 3; nCols = 4;
RspCol = 'b'; CmdCol = 'r--'; mLineWidth = 1.5;

quats = statesArr(7:10,:)';
euls = quat2eul([quats(:, 4), quats(:, 1:3)], 'ZYX'); % convert to scalar first then to euler angles

% Angular velocity
lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,statesArr(11,:)*180/pi,RspCol,tArr,OmegaRefA_BodyArr(1,:)*180/pi,CmdCol, 'LineWidth', mLineWidth);
ylabel('p (deg/s)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,statesArr(12,:)*180/pi,RspCol,tArr,OmegaRefA_BodyArr(2,:)*180/pi,CmdCol, 'LineWidth', mLineWidth);
ylabel('q (deg/s)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,statesArr(13,:)*180/pi,RspCol,tArr,OmegaRefA_BodyArr(3,:)*180/pi,CmdCol, 'LineWidth', mLineWidth);
ylabel('r (deg/s)');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,1),RspCol,tArr,QuatRefA_BodyArr(1,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q1');


lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,2),RspCol,tArr,QuatRefA_BodyArr(2,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q2');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,3),RspCol,tArr,QuatRefA_BodyArr(3,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q3');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
plot(tArr,quats(:,4),RspCol,tArr,QuatRefA_BodyArr(4,:),CmdCol, 'LineWidth', mLineWidth);
ylabel('q4');

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
    leg(r) = plot(tArr,statesArr(13+r,:)*180/pi, 'LineWidth', mLineWidth,'DisplayName',num2str(r));   
end
ylabel('(deg/s)');
xlabel('Time (s)');
legend(leg, 'Interpreter', 'latex');

lm(kk) = subplot(nRows,nCols,kk); kk = kk + 1; hold on;
leg = [];
title('Reaction wheel angular acceleration');
for r = 1:SatellitePlot.params.N_react
    leg(r) = plot(tArr,UarrStore(r,:)*180/pi, 'LineWidth', mLineWidth,'DisplayName',num2str(r));   
end
ylabel('(deg/s2)');
xlabel('Time (s)');
legend(leg, 'Interpreter', 'latex');

linkaxes(lm,'x');