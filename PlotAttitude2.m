
timeDesc = 'Time [s]';
RspCollArr = {[0 0.4470 0.7410],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]	};
LineStyleArr = {'-','-.',':'};
lege = [];
yFontSize = 14;
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
    fprintf('%s Tracking RMS angle error [rad] %f Final enery [J] %f Total compute time [s] %f\n',SatellitePlot.params.Name,RMS_TrackingError_AngRad,Energy2Reac_Total_J_Arr(end),SatellitePlot.TotalRunTime_s);

    if s == 1
        figure(figN); figN = figN + 1; ha = [];
        nRows = 2; nCols = 2;
        RspCol = 'b'; CmdCol = 'r--'; mLineWidth = 1.5;

        tcl = tiledlayout(nRows,nCols,'TileSpacing','tight','Padding','tight'); % "loose", "compact", "tight" or "none"
        % title(tcl,{'MPC tracking',['With basis for GP' num2str(pDimBasisSelect) '''s Y(k-1) and U(k-1) vector elements (subset of full regressor Z(k)']},'FontSize', 10);
    
    end
    RspCol = RspCollArr{s}';
    LineStyle = LineStyleArr{s};

    % Angular velocity
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,w_BA_B_RadsMagArr,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    if s == nSat legeCmd = plot(tArr,w_CmdA_B_RadsMagArr,CmdCol, 'LineWidth', mLineWidth,'DisplayName','Reference'); end
    ylabel('$|w^B|$ [rad/s]','Interpreter','latex','FontSize',yFontSize);

    % angle error
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,Err_BR_AngRadArr,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    ylabel('$\Phi_e$ [rad]','Interpreter','latex','FontSize',yFontSize);

    if false % Plot stats for all wheels
    % Plot reaction wheels
    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    plot(tArr,max(abs(statesArr(13+(1:SatellitePlot.params.N_react),:)),[],1)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);  
    % yline(6500,'r--'); yline(-6500,'r--'); % NanoAvio Max RPM
    ylabel('$|\Omega_i|_{max}$ [RPM]','Interpreter','latex','FontSize',yFontSize);
    xlabel(timeDesc,'FontSize',yFontSize);

    ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    % title('Maximum reaction wheel acceleration');
    lege(end+1) = plot(tArr,max(abs(UarrStore(:,:)),[],1)/(2*pi)*60^2,'Color', RspCol,'linestyle',LineStyle,'LineWidth', mLineWidth,'DisplayName',SatellitePlot.params.Name);  
    MaxTorque = 3.2/1000; % NanoAvio Max torque
    MaxReacWheelAngAccRads2 = MaxTorque/sat_params.I_ws;
    % yline(MaxReacWheelAngAccRads2/(2*pi)*60^2,'r--'); yline(-MaxReacWheelAngAccRads2/(2*pi)*60^2,'r--'); 
    ylabel('$|\dot{\Omega}_i|_{max}$ [RPM$^2$]','Interpreter','latex','FontSize',yFontSize);
    xlabel(timeDesc,'FontSize',yFontSize);

    else
        if true % plot 1 wheel vs the remaining wheels stats
        ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
        plot(tArr,max(abs(statesArr(13+(1:1),:)),[],1)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);  
        % yline(6500,'r--'); yline(-6500,'r--'); % NanoAvio Max RPM
        if s == 1 yline(200,'k--', 'LineWidth', mLineWidth); end;
        ylabel('$|\Omega_1|$ [RPM]','Interpreter','latex','FontSize',yFontSize);
        xlabel(timeDesc,'FontSize',yFontSize);

        ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
        lege(end+1) = plot(tArr,max(abs(statesArr(13+(2:SatellitePlot.params.N_react),:)),[],1)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth,'DisplayName',SatellitePlot.params.Name);  
        % yline(6500,'r--'); yline(-6500,'r--'); % NanoAvio Max RPM
        
        ylabel('$|\Omega_{2,3,4}|_{\mathrm{max}}$ [RPM]','Interpreter','latex','FontSize',yFontSize);
        xlabel(timeDesc,'FontSize',yFontSize);

        else % plot wheel 1 +2 norm vs the remaining wheels stats
        ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
        OmegeReac = statesArr(13+(1:SatellitePlot.params.N_react),:);
        OmegeReac12_Norm = sqrt(OmegeReac(1,:).^2 + OmegeReac(2,:).^2);
        plot(tArr,OmegeReac12_Norm/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);  
        % yline(6500,'r--'); yline(-6500,'r--'); % NanoAvio Max RPM
        if s == 1 yline(200,'k--', 'LineWidth', mLineWidth); end;
        ylabel('$||\Omega_{1,2}||$ [RPM]','Interpreter','latex','FontSize',yFontSize);
        xlabel(timeDesc,'FontSize',yFontSize);

        ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
        lege(end+1) = plot(tArr,max(abs(statesArr(13+(3:SatellitePlot.params.N_react),:)),[],1)/(2*pi)*60,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth,'DisplayName',SatellitePlot.params.Name);  
        % yline(6500,'r--'); yline(-6500,'r--'); % NanoAvio Max RPM        
        ylabel('$|\Omega_{3,4}|_{max}$ [RPM]','Interpreter','latex','FontSize',yFontSize);
        xlabel(timeDesc,'FontSize',yFontSize);    
        end
    end

    % % Energy supplied to reaction wheels
    % ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
    % plot(tArr,Energy2Reac_Total_J_Arr,'Color',RspCol,'linestyle',LineStyle, 'LineWidth', mLineWidth);
    % ylabel('[J]');
    % title('Total energy to reaction wheels','Interpreter','latex');
    
    
end
lg  = legend([lege legeCmd],'Orientation','Horizontal','FontSize',10, 'NumColumns', 3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

linkaxes(ha,'x');