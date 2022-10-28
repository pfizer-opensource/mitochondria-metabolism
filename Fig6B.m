function Fig6B()
load('simpop');

%% Healthy
% Define values of dG and Pi for interpolation
	v_1=sim_pop.pcr_atph(1,:);
	x_1=mean(sim_pop.pih([2:end],:));
% Interpolate values of Pi
	v_2=linspace(-68,-60, 15);  % Delta G range
	pi_healthy=interp1(v_1,x_1,v_2);
% Identify Index at Exhaustion (Pi >= 2)
	ind_h=find(pi_healthy>=2,1);
	pi_healthy(ind_h);
% Interpolate MVO2	
	x_2=(sim_pop.MVO2([2:end],:));
	MVO2_healthy=interp1(v_1,x_2',v_2);
% MVO2 at exhaustion (where Pi >=2)
	Z(1)=mean(MVO2_healthy(ind_h,:));    %Index 
	Z(2)=nan;

%%  Case 1 (DR PDH BetaOxidation)
% Simpop5

clear  v_1 x_1 x_2 v_2
v_2=linspace(-68,-60, 30);  % Delta G range
% Define values of dG and Pi for interpolation
        v_1=sim_pop5.pcr_atph(1,:);
        x_1=mean(sim_pop5.pih([2:end],:));
% Interpolate values of Pi
        pi_case2=interp1(v_1,x_1,v_2);
        ind_c2=find(pi_case2>=2,1);
% Identify Index at Exhaustion (Pi >= 2)
        pi_case2(ind_c2);
% Interpolate MVO2
        x_2=(sim_pop5.MVO2([2:end],:));
        MVO2_c2=interp1(v_1,x_2',v_2);
% MVO2 at exhaustion (where Pi >=2)
        Z(3)=mean(MVO2_c2(ind_c2,:));    

% Add DCA
    v_1=sim_pop6.pcr_atph(1,:);
    x_1=mean(sim_pop6.pih([2:end],:));
% Interpolate values of Pi
    pi_c2dca=interp1(v_1,x_1,v_2);
    ind_c2dca=find(pi_c2dca>=2,1);
% Identify Index at Exhaustion (Pi >= 2)
    pi_c2dca(ind_c2dca);
% Interpolate MVO2
    x_2=(sim_pop6.MVO2([2:end],:));
    MVO2_c2dca=interp1(v_1,x_2',v_2);
% MVO2 at exhaustion (where Pi >=2)
    Z(4)=mean(MVO2_c2dca(ind_c2dca,:)); 



%% Case 2 (Reduceed  mito Volume)  (simpop 7)
clear  v_1 x_1 x_2 v_2
v_2=linspace(-68,-60, 30);  % Delta G range
% Define values of dG and Pi for interpolation
        v_1=sim_pop7.pcr_atph(1,:);
        x_1=mean(sim_pop7.pih([2:end],:));
% Interpolate values of Pi
        pi_c3=interp1(v_1,x_1,v_2);
        ind_c3=find(pi_c3>=2,1);
% Identify Index at Exhaustion (Pi >= 2)
        pi_c3(ind_c3);
% Interpolate MVO2
        x_2=(sim_pop7.MVO2([2:end],:));
        MVO2_c3=interp1(v_1,x_2',v_2);
% MVO2 at exhaustion (where Pi >=2)
        Z(5)=mean(MVO2_c3(ind_c3,:));    

% Add DCA (simpop8)
clear  v_1 x_1 x_2
    v_1=sim_pop8.pcr_atph(1,:);
    x_1=mean(sim_pop8.pih([2:end],:));
% Interpolate values of Pi
    pi_c3dca=interp1(v_1,x_1,v_2);
    ind_c3dca=find(pi_c3dca>=2,1);
% Identify Index at Exhaustion (Pi >= 2)
    pi_c3dca(ind_c3dca);
% Interpolate MVO2
    x_2=(sim_pop8.MVO2([2:end],:));
    MVO2_c3dca=interp1(v_1,x_2',v_2);
% MVO2 at exhaustion (where Pi >=2)
    Z(6)=mean(MVO2_c3dca(ind_c3dca,:)); 
    
%% Case 3 (Decreased nucleotides)
% Simpop 3
clear  v_1 x_1 x_2
% Define values of dG and Pi for interpolation
        v_1=sim_pop3.pcr_atph(1,:);
        x_1=mean(sim_pop3.pih([2:end],:));
% Interpolate values of Pi
        pi_case1=interp1(v_1,x_1,v_2);
        ind_c1=find(pi_case1>=2,1);
% Identify Index at Exhaustion (Pi >= 2)
        pi_case1(ind_c1);
% Interpolate MVO2
        x_2=(sim_pop3.MVO2([2:end],:));
        MVO2_c1=interp1(v_1,x_2',v_2);
% MVO2 at exhaustion (where Pi >=2)
        Z(7)=mean(MVO2_c1(ind_c1,:));    

% Add DCA (Simpop4)
    v_1=sim_pop4.pcr_atph(1,:);
    x_1=mean(sim_pop4.pih([2:end],:));
% Interpolate values of Pi
    pi_c1dca=interp1(v_1,x_1,v_2);
    ind_c1dca=find(pi_c1dca>=2,1);
% Identify Index at Exhaustion (Pi >= 2)
    pi_c1dca(ind_c1dca);
% Interpolate MVO2
    x_2=(sim_pop4.MVO2([2:end],:));
    MVO2_c1dca=interp1(v_1,x_2',v_2);
% MVO2 at exhaustion (where Pi >=2)
    Z(8)=mean(MVO2_c1dca(ind_c1dca,:)); 


%% Plot Figure 
Z2=Z./Z(1).*100;  % Convert to % of healthy MVO2

figure;
colormap(bone)
h=bar(Z2,'grouped','FaceColor','flat'); %ym=bar(z,'grouped','FaceColor',"flat");
hold on;
% Specify face color  
M=bone(8);
inew=1;
for k=1:2:size(Z2,2)   
h.CData(k,:) =M(inew,:); % Healthy
inew=inew+1;
end
inew=4;
for k=2:2:size(Z2,2)
h.CData(k,:) =M(inew,:); % Healthy
inew=inew+1;
end
sn= length(MVO2_healthy(ind_h,:));
std2=zeros(sn,1);
vv=[MVO2_healthy(ind_h,:)',std2, MVO2_c2(ind_c2,:)', MVO2_c2dca(ind_c2dca,:)', MVO2_c3(ind_c3,:)',...
MVO2_c3dca(ind_c3dca,:)', MVO2_c1(ind_c1,:)',MVO2_c1dca(ind_c1dca,:)'];
clear ub lb
for kk=1:length(vv(1,:))
    ub(kk)=prctile(vv(:,kk),95)-(Z(kk));
    lb(kk)=-prctile(vv(:,kk),5)+(Z(kk));
end
errorbar([1:8], Z2, lb.*100, ub.*100, 'k', 'linestyle', 'none');


ylim([0  109]); ylabel('VO_2 (% of Healthy)')
title('Simulated Oxygen Consumption at Fatigue')
xtickangle(45);
ax=gca;
ax.XTickLabel={'Healthy','','\downarrow  PDH & \beta-Ox'...
'\downarrow  PDH & \beta-Ox + DCA','\downarrow Mito Volume',...
'\downarrow Mito Volulme + DCA','\downarrow Nuceotides','\downarrow Nucleotides +  DCA'};

saveas(gcf,'Fig6B.pdf');
