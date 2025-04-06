
timeDesc = 'Time [s]';
RspCollArr = {[0 0.4470 0.7410],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]	};
lege = [];
yFontSize = 12;
for s = 1:nSat
     kk = 1;
    
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

    omega_reac = SatellitePlot.statesArr(13+(1:SatellitePlot.params.N_react),2:end); %statesArr(,:);
    Power_reac = SatellitePlot.params.I_ws*omega_reac.*UarrStore;
    Energy2Reac_Total_J_Arr = cumtrapz(sum(Power_reac*Ts,1));

    RMS_TrackingError_AngRad = sqrt(sum(Err_BR_AngRadArr.^2)/Nsim);
    fprintf('%s Tracking RMS angle error [deg] %f Final enery [J] %f Total compute time [s] %f\n',SatellitePlot.params.Name,RMS_TrackingError_AngRad*180/pi,Energy2Reac_Total_J_Arr(end),SatellitePlot.TotalRunTime_s);

    if s == 1
        figure(figN); figN = figN + 1; ha = [];
        nRows = 4; nCols = 2;
        RspCol = 'b'; CmdCol = 'r--'; mLineWidth = 1.5;

        tcl = tiledlayout(nRows,nCols,'TileSpacing','tight','Padding','tight'); % "loose", "compact", "tight" or "none"
        % title(tcl,{'MPC tracking',['With basis for GP' num2str(pDimBasisSelect) '''s Y(k-1) and U(k-1) vector elements (subset of full regressor Z(k)']},'FontSize', 10);
    
    end
    RspCol = RspCollArr{s}';
    LineStyle = LineStyleArr{s};
    
    % quat
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,quats(:,1),'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,QuatRefA_A_Arr(1,:),CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$q_1$','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,quats(:,2),'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,QuatRefA_A_Arr(2,:),CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$q_2$','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,quats(:,3),'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,QuatRefA_A_Arr(3,:),CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$q_3$','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,quats(:,4),'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,QuatRefA_A_Arr(4,:),CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$q_4$','Interpreter','latex','FontSize',yFontSize);

    % Angular velocity
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,pRadsArr,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,OmegaCmdA_BodyArr(1,:),CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$w_1$ [rad/s]','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,qRadsArr,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,OmegaCmdA_BodyArr(2,:),CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$w_2$ [rad/s]','Interpreter','latex','FontSize',yFontSize);


    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,rRadsArr,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,OmegaCmdA_BodyArr(3,:),CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$w_3$ [rad/s]','Interpreter','latex','FontSize',yFontSize);
    xlabel(timeDesc,'FontSize',yFontSize);


    % Angular velocity
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    lege(end+1) = plot(tArr,w_BA_B_RadsMagArr,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth,'DisplayName',SatellitePlot.params.Name);
    if s == nSat legeCmd = plot(tArr,w_CmdA_B_RadsMagArr,CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$|w^B|$ [rad/s]','Interpreter','latex','FontSize',yFontSize);
    xlabel(timeDesc,'FontSize',yFontSize);

    
end
lg  = legend([lege legeCmd],'Orientation','Horizontal','FontSize',10, 'NumColumns', 3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

linkaxes(ha,'x');

