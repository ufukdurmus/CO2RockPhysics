% Function to make well logs plotting easier and prettier
% coded by Ufuk Durmus on 7/7/2023
% should be in the same directory
% variable names should be in the correct order
% example ylimit -> %ylim -> depth axis -> ylimit = [300 1380]; -> ylimit=[-inf inf];

function well_logs_disp_final(depth,depth2,depth_core,ylimit,cali,gr,vp,vs,rhob,rhob_core,drho,phit,phit_core,vclay,intra_drake,lower_drake1,cook4,cook2,burton,johansen4,johansen2,amundsen)

%% CALI
figure
subplot(1,8,1)
plot(cali,depth,'r','LineWidth',2)
set(gca,'YDir','reverse')
xlim([6 12])
ylim(ylimit)
xlabel('CALI (in)','FontSize',12,'FontWeight','bold')
ylabel('Depth (m)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';

%% GR
subplot(1,8,2)
plot(gr,depth,'g','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlim([20 200])
ylim(ylimit)
xlabel('GR (api)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';

%% Vp
subplot(1,8,3)
plot(vp,depth,'k','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
ylim(ylimit)
xlim([1000 inf])
xlabel('P-wave (m/s)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';

%% Vs
subplot(1,8,4)
plot(vs,depth,'b','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
ylim(ylimit)
xlim([500 inf])
xlabel('S-wave (m/s)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';

%% RHOB
subplot(1,8,5)
scatter(rhob_core,depth_core,15,'k','filled')
hold on
plot(rhob,depth,'LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
ylim(ylimit)
xlim([1.85 2.85])
xlabel('RHOB (g/cc)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.ColorOrder = [0.8500 0.3250 0.0980];
ax.XGrid = 'on';
ax.YGrid = 'off';

%% DRHO
subplot(1,8,6)
plot(drho,depth,'LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
ylim(ylimit)
xlim([-0.1 0.1])
xlabel('DRHO','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.ColorOrder = [0.4940 0.1840 0.5560];
ax.XGrid = 'on';
ax.YGrid = 'off';

%% PHIT
% Clasic Display
subplot(1,8,7)
scatter(phit_core,depth_core,15,'k','filled')
hold on
plot(phit,depth2,'LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
ylim(ylimit)
xlim([0 0.4])
xlabel('Porosity (v/v)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.ColorOrder = [0.9290 0.6940 0.1250];
ax.XGrid = 'on';
ax.YGrid = 'off';

% visualize the tops
r_tops = 0:0.1:1; % showing the tops in the range as in vclay section
%% Volumetric Mineral
subplot(1,8,8)
plot(vclay,depth2,'LineWidth',2)
hold on
plot(r_tops,intra_drake.*ones(size(r_tops)),'black','LineWidth',2)
text(1,intra_drake,'intra-drake','Color','black','FontWeight','bold')
hold on
plot(r_tops,lower_drake1.*ones(size(r_tops)),'black','LineWidth',2)
text(1,lower_drake1,'lower-drake1','Color','black','FontWeight','bold')
hold on
plot(r_tops,cook4.*ones(size(r_tops)),'black','LineWidth',2)
text(1,cook4-2,'cook4','Color','black','FontWeight','bold')
hold on
plot(r_tops,cook2.*ones(size(r_tops)),'black','LineWidth',2)
text(1,cook2+1,'cook2','Color','black','FontWeight','bold')
hold on
plot(r_tops,burton.*ones(size(r_tops)),'black','LineWidth',2)
text(1,burton,'burton','Color','black','FontWeight','bold')
hold on
plot(r_tops,johansen4.*ones(size(r_tops)),'black','LineWidth',2)
text(1,johansen4,'johansen4','Color','black','FontWeight','bold')
hold on
plot(r_tops,johansen2.*ones(size(r_tops)),'black','LineWidth',2)
text(1,johansen2,'johansen2','Color','black','FontWeight','bold')
hold on
plot(r_tops,amundsen.*ones(size(r_tops)),'black','LineWidth',2)
text(1,amundsen,'amundsen','Color','black','FontWeight','bold')
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
ylim(ylimit)
xlim([0 1])
xlabel('Vclay (v/v)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.ColorOrder = [0.6350 0.0780 0.1840];
ax.XGrid = 'on';
ax.YGrid = 'off';

end
