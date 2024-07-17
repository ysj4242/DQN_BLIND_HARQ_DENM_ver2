close all    % Close all open figures
clear        % Reset variables
clc          % Clear the command window

T =50;

time = (1e-3:1e-3:T);  %%% 0.001간격
time = reshape(time,length(time),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 300
A = load('6.mat');        %black
B = load('7.mat');        %
C = load('8.mat');         %red
D = load('9.mat');         %
E = load('8.mat');         %blue
F = load('8.mat');         %
G = load('8.mat');         %green
H = load('8.mat');         %
I = load('8.mat');          %green

% C = load('control_524_V2I_I2V_delt.mat');
 
% A = load('2rep_800_np_d1_0.8.mat');         %black
% B = load('2rep_800_np_d1_0.65.mat');        %black
% C = load('2rep_800_np_nd.mat');             %black
% D = load('2rep_400_np_d1_0.65.mat');        %red
% E = load('2rep_400_np_d1_0.3.mat');         %red
% F = load('2rep_400_np_nd.mat');             %red
% G = load('2rep_200_np_d1_0.65.mat');        %green
% H = load('2rep_200_np_d1_0.3.mat');         %green
% I = load('2rep_200_np_nd.mat');         %green

% A = load('2rep_800_np_pd1.mat');         %black
% B = load('2rep_800_np_d1_1.mat');         %red
% C = load('2rep_800_np_nd.mat');           %blue
% D = load('2rep_200_np_pd1.mat');          %green
% E = load('2rep_200_np_d1_1.mat');         %magenta
% F = load('2rep_200_p5_dc.mat');           %
% G = load('2rep_200_p5_nd.mat');           %
% H = load('2rep_200_np_nd.mat');           %

% a1 = A.outputValues.packetReceptionRatioCV2X;
% b1 = B.outputValues.packetReceptionRatioCV2X;
% c1 = C.outputValues.packetReceptionRatioCV2X;
% d1 = D.outputValues.packetReceptionRatioCV2X;
% e1 = E.outputValues.packetReceptionRatioCV2X;
% f1 = F.outputValues.packetReceptionRatioCV2X;
% g1 = G.outputValues.packetReceptionRatioCV2X;
% h1 = H.outputValues.packetReceptionRatioCV2X;
% i1 = H.outputValues.packetReceptionRatioCV2X;

% 
a1 = A.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
b1 = B.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
c1 = C.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
d1 = D.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
e1 = E.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
f1 = F.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
g1 = G.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
h1 = H.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);
i1 = H.outputValues.packetReceptionRatioCV2X_DENM(:,:,[1,3,5,6,7]);

% a1 = A.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% b1 = B.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% c1 = C.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% d1 = D.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% e1 = E.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% f1 = F.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% g1 = G.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% h1 = H.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);
% i1 = H.outputValues.packetReceptionRatioCV2X_CAM(:,:,[1,3,5,6,7]);

a2 = reshape(a1,1,[]);
b2 = reshape(b1,1,[]);
c2 = reshape(c1,1,[]);
d2 = reshape(d1,1,[]);
e2 = reshape(e1,1,[]);
f2 = reshape(f1,1,[]);
g2 = reshape(g1,1,[]);
h2 = reshape(h1,1,[]);
i2 = reshape(i1,1,[]);

rep2_a = [1,a2];
rep1_b = [1,b2];
rep2_c = [1,c2];
rep1_d = [1,d2];
rep2_e = [1,e2];
rep1_f = [1,f2];
rep2_g = [1,g2];
rep1_h = [1,h2];
rep2_i = [1,i2];

% LDM_DCC_bak = (Ncontrol - control_cor)./((Ncontrol+control_cor)/2).*100;
% LDM_DCC_NR_V2X_200 = [1,c2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 500
% E = load('NO_DCC_NR_500.mat');
% F = load('ETSI_DCC_NR_500.mat');
% G = load('LDM_DCC_NR_500.mat');
% 
% a1 = E.outputValues.packetReceptionRatioCV2X;
% b1 = F.outputValues.packetReceptionRatioCV2X;
% c1 = G.outputValues.packetReceptionRatioCV2X;
% 
% e2 = reshape(a1,1,10);
% f2 = reshape(b1,1,10);
% g2 = reshape(c1,1,10);
% 
% NO_DCC_NR_V2X_500 = [1,e2];
% DCC_ACTIVE_NR_V2X_500 = [1,f2];
% LDM_DCC_NR_V2X_500 = [1,g2];
% 
seung = [50 100 150 200 300];
distance = [0, seung];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
semilogy(distance,rep2_a,'Color','black','LineStyle','-','Marker','o','LineWidth',1.5);
hold on
semilogy(distance,rep1_b,'Color','black','LineStyle','--','Marker','o','LineWidth',1.5);
hold on
semilogy(distance,rep2_c,'Color','r','LineStyle','-','Marker','o','LineWidth',1.5);
hold on
semilogy(distance,rep1_d,'Color','r','LineStyle','--','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep2_e,'Color','blue','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_f,'Color','blue','LineStyle','--','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep2_g,'Color','green','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_h,'Color','green','LineStyle','--','Marker','o','LineWidth',1.5);
% semilogy(distance,DCC_ACTIVE_NR_V2X_500,'Color','r','LineStyle','-','Marker','x','LineWidth',1.5);
% semilogy(distance,LDM_DCC_NR_V2X_500,'Color','r','LineStyle','-','Marker','square','LineWidth',1.5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% semilogy(distance,rep2_a,'Color','black','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_b,'Color','black','LineStyle',':','Marker','diamond','LineWidth',1.5);
% hold on
% semilogy(distance,rep2_c,'Color','black','LineStyle','--','Marker','*','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_d,'Color','red','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep2_e,'Color','red','LineStyle',':','Marker','diamond','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_f,'Color','red','LineStyle','--','Marker','*','LineWidth',1.5);
% hold on
% semilogy(distance,rep2_g,'Color','green','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_h,'Color','green','LineStyle',':','Marker','diamond','LineWidth',1.5);
% hold on
% semilogy(distance,rep2_i,'Color','green','LineStyle','--','Marker','*','LineWidth',1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% semilogy(distance,rep2_a,'Color','black','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_b,'Color','red','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep2_c,'Color','blue','LineStyle','-','Marker','o','LineWidth',1.5);
% hold on
% semilogy(distance,rep1_d,'Color','green','LineStyle','-','Marker','o','LineWidth',3);
% hold on
% semilogy(distance,rep2_e,'Color','m','LineStyle','-','Marker','o','LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


grid on
% lgd = legend({'rep=2, ρ=800','rep=1, ρ=800','rep=2, ρ=400','rep=1, ρ=400','rep=2, ρ=200','rep=1, ρ=200'},'Orientation','horizontal');
% lgd = legend({'p5, ρ=800','np, ρ=800','p5, ρ=400','np, ρ=400','p5, ρ=200','np, ρ=200',},'Orientation','horizontal');
% lgd = legend({'p20, ρ=400','np, ρ=400','p20, ρ=200','np, ρ=200'},'Orientation','horizontal');
% lgd = legend({'ρ=400, df=0.01','ρ=400, df=0.1','ρ=400, df=0.5','ρ=400, df=1','ρ=400, df=1.5'},'Orientation','horizontal');

% lgd.NumColumns = 1;

% lgd = legend({'pd, ρ=800','d1, ρ=800','nd, ρ=800','pd, ρ=400','d1, ρ=400','nd, ρ=400','pd, ρ=200','d1, ρ=200','nd, ρ=200'},'Orientation','horizontal');
lgd = legend({'adaptive DENM, ρ=800','fixed DENM, ρ=800','adaptive DENM, ρ=400','fixed DENM, ρ=400','adaptive DENM, ρ=200','fixed DENM, ρ=200'},'Orientation','horizontal');
lgd.NumColumns = 2;

xlabel('Distance');
ylabel('Packet Delivery Ratio (PDR)');
% title('PDR');
title('DENM PDR');
% title('CAM PDR');
axis([0 300 0 1])
set(gca,'FontSize',12)
set(gca,'YTick', [0:0.1:1])

% figure(2)
% semilogy(time,A.sinrManagement.meanCBR,'linewidth',1,'Color','black');
% hold on
% semilogy(time,B.sinrManagement.meanCBR,'linewidth',1,'Color','r');
% grid on
% legend('Route','non-route');
% xlabel('Simulation time (ms)','FontSize',12,'FontWeight','bold','Color','k')
% ylabel('Mean CBR','FontSize',12,'FontWeight','bold','Color','k');
% title('Computation for all VUEs');
% set(gca,'FontSize',12)
% axis([0 T 0.15 1])