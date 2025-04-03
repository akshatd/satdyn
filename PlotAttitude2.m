

RspCollArr = {[0 0.4470 0.7410],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]	};
lege = [];
for s = 1:nSat
     kk = 1;
    if s == 1
        figure(figN); figN = figN + 1; ha = [];
        nRows = 2; nCols = 2;
        RspCol = 'b'; CmdCol = 'r--'; mLineWidth = 1.5;

        tcl = tiledlayout(nRows,nCols,'TileSpacing','tight','Padding','tight'); % "loose", "compact", "tight" or "none"
        % title(tcl,{'MPC tracking',['With basis for GP' num2str(pDimBasisSelect) '''s Y(k-1) and U(k-1) vector elements (subset of full regressor Z(k)']},'FontSize', 10);
    
    end
    RspCol = RspCollArr{s}';

    SatellitePlot = SatelliteArr(s); % select satellite to plot
    statesArr = SatellitePlot.statesArr(:,1:end-1);
    UarrStore = SatellitePlot.UarrStore;
    tArr = SatellitePlot.tArr;
    OmegaCmdA_BodyArr = SatellitePlot.OmegaCmdA_BodyArr;
    QuatRefA_A_Arr = SatellitePlot.QuatRefA_A_Arr;
    Err_BR_AngRadArr = SatellitePlot.Err_BR_AngRadArr;

    quats = statesArr(7:10,:)';
    euls = quat2eul([quats(:, 4), quats(:, 1:3)], 'ZYX'); % convert to scalar first then to euler angles
    
    pRadsArr= statesArr(11,:); qRadsArr= statesArr(12,:); rRadsArr= statesArr(13,:);
    w_BA_B_RadsMagArr = sqrt(pRadsArr.^2 + qRadsArr.^2 + rRadsArr.^2);
    w_CmdA_B_RadsMagArr = sqrt(OmegaCmdA_BodyArr(1,:).^2 + OmegaCmdA_BodyArr(2,:).^2 + OmegaCmdA_BodyArr(3,:).^2);

    % Angular velocity
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,w_BA_B_RadsMagArr*180/pi,'Color',RspCol, 'LineWidth', mLineWidth);
    if s == 1 legeCmd = plot(tArr,w_CmdA_B_RadsMagArr*180/pi,CmdCol, 'LineWidth', mLineWidth,'DisplayName','Command'); end
    ylabel('(deg/s)');
    title('Angular velocity $w_{BA}^B$ mag','Interpreter','latex');

    % angle error
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,Err_BR_AngRadArr*180/pi,'Color',RspCol, 'LineWidth', mLineWidth);
    ylabel('[deg]');
    title('Angle error ','Interpreter','latex');

    % Plot reaction wheels
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,max(abs(statesArr(13+(1:SatellitePlot.params.N_react),:)),[],1)/(2*pi)*60,'Color',RspCol, 'LineWidth', mLineWidth);  
    % yline(6500,'r--'); yline(-6500,'r--'); % NanoAvio Max RPM
    ylabel('(RPM)');
    xlabel('Time (s)');
    title('Maximum reaction wheel speed','Interpreter','latex');

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    title('Maximum reaction wheel acceleration');
    lege(end+1) = plot(tArr,max(abs(UarrStore(:,:)),[],1)/(2*pi)*60^2,'Color', RspCol,'LineWidth', mLineWidth,'DisplayName',SatellitePlot.params.Name);  
    MaxTorque = 3.2/1000; % NanoAvio Max torque
    MaxReacWheelAngAccRads2 = MaxTorque/sat_params.I_ws;
    % yline(MaxReacWheelAngAccRads2/(2*pi)*60^2,'r--'); yline(-MaxReacWheelAngAccRads2/(2*pi)*60^2,'r--'); 
    ylabel('(RPM2)');
    xlabel('Time (s)');
    
end
lg  = legend([lege legeCmd],'Orientation','Horizontal','FontSize',10, 'NumColumns', 2); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout
