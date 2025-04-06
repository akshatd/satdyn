
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

    % Reaction wheel speed
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,omega_reac(1,:)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\Omega_1$ [RPM]','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,omega_reac(2,:)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\Omega_2$ [RPM]','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,omega_reac(3,:)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\Omega_3$ [RPM]','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,omega_reac(4,:)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\Omega_4$ [RPM]','Interpreter','latex','FontSize',yFontSize);

    % Reaction wheel acceleration
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,UarrStore(1,:)/(2*pi)*60^2,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\dot{\Omega}_1$ [RPM$^2$]','Interpreter','latex','FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,UarrStore(2,:)/(2*pi)*60^2,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\dot{\Omega}_2$ [RPM$^2$]','Interpreter','latex','FontSize',yFontSize);


    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,UarrStore(3,:)/(2*pi)*60^2,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\dot{\Omega}_3$ [RPM$^2$]','Interpreter','latex','FontSize',yFontSize);
    xlabel(timeDesc,'FontSize',yFontSize);
    
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    lege(end+1) = plot(tArr,UarrStore(4,:)/(2*pi)*60^2,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth,'DisplayName',SatellitePlot.params.Name);
    ylabel('$\dot{\Omega}_4$ [RPM$^2$]','Interpreter','latex','FontSize',yFontSize);
    xlabel(timeDesc,'FontSize',yFontSize);
end
lg  = legend([lege],'Orientation','Horizontal','FontSize',10, 'NumColumns', 3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

linkaxes(ha,'x');

