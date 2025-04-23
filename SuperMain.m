clear all
Omega1_LimitArr_RPM = 0:50:450;
RMS_TrackingError_AngRad_Arr = [];
Omega234_Max_Max_Arr = [];
for li = 1:length(Omega1_LimitArr_RPM)
    Omega1_Limit_RPM = Omega1_LimitArr_RPM(li);

    main;

    RMS_TrackingError_AngRad_Arr(li) = RMS_TrackingError_AngRad;

    Omega234_Max_Max_Arr(li) = max(max(abs(statesArr(13+(2:SatellitePlot.params.N_react),:)),[],1)/(2*pi)*60); % max omega2,3,4 over all time
end

%%
figure(20);  ha = []; kk =1;
nRows = 1; nCols = 2;
mLineWidth = 1.5;
RspCollArr = {[0 0.4470 0.7410],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]	};
yFontSize = 14;
tcl = tiledlayout(nRows,nCols,'TileSpacing','tight','Padding','tight'); % "loose", "compact", "tight" or "none"

ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
plot(Omega1_LimitArr_RPM,RMS_TrackingError_AngRad_Arr,'Color',RspCollArr{1}','LineWidth',1.5); 
plot(450,0.0191,'o','Color',RspCollArr{1}','LineWidth',1.5)
yline(0.0204,'r--','LineWidth',1.5)
ylabel('RMS $\Phi_e$ [rad]','Interpreter','latex','FontSize',yFontSize)
xlabel('$|\Omega_1|$ limit [RPM]','Interpreter','latex','FontSize',yFontSize)
legend('NMPCdeg','NMPCnom','PDnom and PDdeg','Interpreter','latex','FontSize',10)
xlim([0 450])
% textbox NMPCnom peak* |\Omega_1|

ha(kk) = nexttile(kk); kk = kk + 1; hold on; 
plot(Omega1_LimitArr_RPM,Omega234_Max_Max_Arr,'Color',RspCollArr{1}','LineWidth',1.5); 
plot(445,434.5,'o','Color',RspCollArr{1}','LineWidth',1.5)
yline(442,'r--','LineWidth',1.5)
yline(834,'k--','LineWidth',1.5)
ylabel('$\mathrm{Peak}^{*} \mathrm{  } \mathrm{  } |\Omega_{2,3,4}|_{\mathrm{max}}$ [RPM]','Interpreter','latex','FontSize',yFontSize)
xlabel('$|\Omega_1|$ limit [RPM]','Interpreter','latex','FontSize',yFontSize)
legend('NMPCdeg','NMPCnom','PDnom','PDdeg','Interpreter','latex','FontSize',10)
xlim([0 450])