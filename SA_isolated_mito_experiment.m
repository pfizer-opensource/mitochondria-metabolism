% This is the main driver program for simulating steady-state cardiac
% energetics of the control hearts under different work rates.
% Written by Fan Wu and Daniel Beard, Medical College of Wisconsin, 2008

setup1;
param(39) = 250; % leak parameter
param(1)=param(1)/100; % LFM Changed 2/11/2020 for pdh model substrate selection    PDH . - 2.05e-5 M/s
param(2)= param(2)/100;% LFM Changed 2/11/2020 for pdh model substrate selction   CS . - 9.82e-4 M/s
param(43) = 25e-6; %4e-4; %Carnitine in the Cytosol (400 uM)
param(44)= 5.6e-5; %ki of CnAcyl carnitines in the cytosol for CACT 56 uM
param(45)= 2.0e-4; %ki of CACT for carnitine in the cytosol


   load('sf_out')
   load('Varstruc_default1.mat')
   delta_sf=sf_out;

global iPO2 idPsi iH_x iATP_x iADP_x iAMP_x iGTP_x iGDP_x iPI_x iNADH_x iQH2_x ...
       iOAA_x iACCOA_x iCIT_x iICIT_x iAKG_x iSCOA_x iCOASH_x iSUC_x iFUM_x ...
       iMAL_x iGLU_x iASP_x iK_x iMg_x iCO2tot_x iCred_i iATP_i iADP_i iAMP_i ...
       iPI_i iH_i iMg_i iK_i iATP_c iADP_c iPI_c iH_c iMg_c iK_c iPYR_x iPYR_i ...
       iPYR_c iCIT_i iCIT_c iAKG_i iAKG_c iSUC_i iSUC_c iMAL_i iMAL_c iASP_i ...
       iASP_c iGLU_i iGLU_c iFUM_i iFUM_c iICIT_i iICIT_c iPCr_c iAMP_c iCr_c;
     
global i_J_C1 i_J_C3 i_J_C4 i_J_F1 i_J_ANT i_J_PI1 i_J_Hle i_J_KH i_J_pdh i_J_cits ...
       i_J_acon i_J_isod i_J_akgd i_J_scoas i_J_sdh i_J_fum i_J_mdh i_J_ndk i_J_got ...
	   i_J_PI2 i_J_ATPt i_J_ADPt i_J_AMPt i_J_PYR_H i_J_GLU_H i_J_CIT_MAL i_J_SUC_MAL ... 
       i_J_AKG_MAL i_J_MAL_PI i_J_ASP_GLU i_J_PYRt i_J_CITt i_J_ICITt i_J_MALt ... 
	   i_J_AKGt i_J_SUCt i_J_FUMt i_J_GLUt i_J_ASPt i_J_CKe i_J_AKe i_J_AtC...
       iC16Carn_m iC16CoA_m iC16EnoylCoA_m iC16OHCoA_m iC16KetoCoA_m...   
       iC14Carn_cy iC14Carn_m iC14CoA_m iC14EnoylCoA_m iC14OHCoA_m...    
       iC14KetoCoA_m iC12Carn_cy iC12Carn_m iC12CoA_m iC12EnoylCoA_m iC12OHCoA_m...   
       iC12KetoCoA_m iC10Carn_cy iC10Carn_m iC10CoA_m iC10EnoylCoA_m iC10OHCoA_m...  
       iC10KetoCoA_m iC8Carn_cy iC8Carn_m iC8CoA_m iC8EnoylCoA_m iC8OHCoA_m...  
       iC8KetoCoA_m iC6Carn_cy iC6Carn_m iC6CoA_m iC6EnoylCoA_m iC6OHCoA_m...
       iC6KetoCoA_m iC4Carn_cy iC4Carn_m iC4CoA_m...
       iC4EnoylCoA_m iC4OHCoA_m iC4KetoCoA_m iAcetylCoAMAT iFADH_m iNADHm iCoAMAT...
       iC16AcylCoACYT iC16Carn_cy;

options = odeset('NonNegative', 3:109);
xo=xo_con;
xo(63:109)=0;
%%
%   Pyr_c;  C16carn_c; C14carn_c; C12; C10;   C8;    C6;     C4;
IC=[7.5e-5; 1.71e-7; 2.3e-8; 1.1e-7; 1.9e-8; 5.2e-8; 1.7e-8; 8e-9; 2.68e-5;];  


%%
SF=10^6; %M
xo(105) = 0; 
xo(106) = 0.46/SF;  % FADH2 Matrix x(44)
xo(108) = 0;% 4979.0/SF; %CoASH Rule x(46)
xo(110) = 25/SF; % 25 uM Carnitine in cytosol
xo(111) = 25/SF; %25 uM Carnitine in the intermembrane space

% indices:[cact_f; cact_r(C4-C16);cpt2 (C4-C16); Vvlcad (C12-16); Vlcad(c8-C16);  Vmcad (C4-C6);   Vscad (C4-C6); Vcrot; Vmschad;  Vmckat; Vmpt; Vcpt1 (C16); ]
%         [ 1      2;               3;             4;                 5;              6;              7;            8;      9;      10;     11;         12;]  
% Original Beta oxidation Vmax parameters

     sig=0.5;                                         % Variability
     nm = 20;                                         % # of Parameter Sets in LHS
     Ind=[5; 9; 23; 35; 39; 12; 1; ];                 % Indices for parameters (VmaxAKG, VmaxMDH, Ki_NADHH_akg, ANT activity,Proton Leak Activity, Pyruvate Transporter, VmaxPDH, CPT1 transporter, KmCpt1) 
     nInd = numel(Ind);
     sf_out=sf_out';
     sfInd=numel(sf_out); %sf_out + sf_out.*sig %param(Ind(1:nInd)) - param(Ind(1:nInd)).*sig
     
     ub= [param(Ind(1:nInd)) + param(Ind(1:nInd)).*sig, Varstruc.boxVmax(12)+ Varstruc.boxVmax(12).*sig, Varstruc.Kmcpt1CarCYT + Varstruc.Kmcpt1CarCYT.*sig ];%        [param([Ind(1:nInd)]) + param([Ind(1:nInd)]).*sig];%, 
     lb= [param(Ind(1:nInd)) - param(Ind(1:nInd)).*sig, Varstruc.boxVmax(12)- Varstruc.boxVmax(12).*sig, Varstruc.Kmcpt1CarCYT - Varstruc.Kmcpt1CarCYT.*sig ]; %[param([Ind(1:nInd)]) - param([Ind(1:nInd)]).*sig];%, 
     X_lhs=lhsdesign(nm,nInd+2);
     X_2= bsxfun(@plus,lb,bsxfun(@times,X_lhs,(ub-lb)));


nInd=numel(Ind);
sim_JO2_Carn=ones(length(X_2(:,1)),7);

tic
for k=1:1:length(X_2(:,1))
        
        param1=param;
        param1(Ind(1:7)) =X_2(k,1:7); %norm_guess(param(m),param(m)*0.2);
        Varstruc.boxVmax(12)=X_2(k,8); 
        Varstruc.Kmcpt1CarCYT =X_2(k,9);  
        boxVmax=Varstruc.boxVmax;

%% Set up SA
 % Parmvary= [ 0.01; 0.1; 1; 5; 10; 100;];
 

     


%% Set up the initial conditions
% Run the simulation to an initial state 

xo(iACCOA_x) = 70e-6; 
xo(iNADH_x) = 16/SF; 
xo([42:62 iATP_c iADP_c iAMP_c iATP_i iADP_i iAMP_i]) = 0;
[t0,x0] = ode15s(@Mito_dXdT,[0 60],xo,options,param1,boxVmax,delta_sf);

% Do state-1 experiment with experimental buffer levels
xo1 = x0(end,:);
xo1(iPI_c) = 5.0e-3;
xo1(iH_c) = 10.^(-7.2);
[t1,x1] = ode15s(@Mito_dXdT,[0 60],xo1,options,param1,boxVmax,delta_sf);

    Carnind = [63,69,75,81,87,93,99]; %63
    J_O2_el=ones(numel(Carnind),1);
    for i = 1:numel(Carnind)
    % Add substrates and run state-2 simulation
    xo2 = x1(end,:);
    xo2(iMAL_c) = 1.0e-3;
    xo2(Carnind(i))=20e-6; 
    [t2,x2] = ode15s(@Mito_dXdT,[0 60],xo2,options,param1,boxVmax,delta_sf);

    % Add ADP for state 3 simulation
    xo3 = x2(end,:);
    xo3(iADP_c) = 0.375e-3;
    [t3,x3] = ode15s(@Mito_dXdT,[61 1200],xo3,options,param1,boxVmax,delta_sf);

    clear t
    clear x
    % Combine simulation timecourses
    t = [t2(2:end); t3(2:end)];
    x = [x2(2:end,:); x3(2:end,:)];

    clear Flux
    clear J_O2
        for j = 1:length(t)
        [Flux, Fluxbox(j,:)] = Mito_Flux(x(j,:),param1,boxVmax,delta_sf);
        J_O2(j) = Flux(6)/2; % J_O2 ( mol / sec / L mito )         
        end

        clear J_electrode

% Accounting for electrode response time
    [t_O2,J_electrode] = ode15s(@dxdt_electrode, t, 0, [], t, J_O2);  
    tss=find(t>90,1);
 
    J_O2_el(i) = J_electrode(tss); % max(J_electrode);
     
    end
    
    sim_JO2_Carn(k,:)=J_O2_el;
end
    toc


        
 
%% Parameters varied to explore parameter space and expected variability
% Vmax of the beta-oxidation parameters is varied by +/- 10% 
param(43)=0;
% for m= 1:20 
% boxVmax1 = [7e-9; 7e-9; 6.5167e-9; 1.3333e-10; 1.6667e-10;  1.35e-9;  1.35e-9;  6.0e-8;  1.6667e-8;  6.2833e-9;  4.7333e-8; 2e-10;]; % Bakker Vmax values (M/s)
% boxmin=boxVmax1 - (boxVmax1*0.1);
% boxmax=boxVmax1 + (boxVmax1*0.1);
% boxVmax=boxmin + rand(1,12)'.*(boxmax-boxmin);
sim_JO2_nocarn=ones(length(X_2(:,1)),7);
for kk=1:1:length(X_2(:,1))
 
        
        param1=param;
        param1(Ind(1:7)) =X_2(kk,1:7); 
        %sf_out=(X_2(k,1:sfInd));

        Varstruc.boxVmax(12)=X_2(kk,8);  %n+5%Varstruc.boxVmax(n)=norm_guess(Varstruc.boxVmax(n),Varstruc.boxVmax(n)*0.1)
        Varstruc.Kmcpt1CarCYT =X_2(kk,9);%nn+6
        boxVmax=Varstruc.boxVmax;


%% Set up SA
 % Parmvary= [ 0.01; 0.1; 1; 5; 10; 100;];
  %ind=[44,45];%  sf indices
 %for h= 1: numel(Parmvary);
      %param1=param;
  %   param1(44)=(param(44)*Parmvary(h));

%% Set up the initial conditions
% Run the simulation to an initial state 
%xo(iCOASH_x) = 5E-3;
xo(iACCOA_x) = 70e-6; 
xo(iNADH_x) = 16/SF; 
xo([42:62 iATP_c iADP_c iAMP_c iATP_i iADP_i iAMP_i]) = 0;
[t0,x0] = ode15s(@Mito_dXdT,[0 60],xo,options,param1,boxVmax,delta_sf);

% Do state-1 experiment with experimental buffer levels
xo1 = x0(end,:);
xo1(iPI_c) = 5.0e-3;
xo1(iH_c) = 10.^(-7.2);
[t1,x1] = ode15s(@Mito_dXdT,[0 60],xo1,options,param1,boxVmax,delta_sf);


  %  Parmvary=[0.01; 0.1; 1; 10; 100;];
    Carnind = [63,69,75,81,87,93,99]; %63
    J_O2_el_n=ones(numel(Carnind),1);
    for ii = 1:numel(Carnind)

% Add substrates and run state-2 simulation
    xo2 = x1(end,:);
    xo2(iMAL_c) = 1.0e-3;
    xo2(Carnind(ii))=20e-6; 
    [t2,x2] = ode15s(@Mito_dXdT,[0 60],xo2,options,param1,boxVmax,delta_sf);

% Add ADP for state 3 simulation
    xo3 = x2(end,:);
    xo3(iADP_c) = 0.375e-3;
    [t3,x3] = ode15s(@Mito_dXdT,[61 1200],xo3,options,param1,boxVmax,delta_sf);

    clear t
    clear x
% Combine simulation timecourses
    t = [t2(2:end); t3(2:end)];
    x = [x2(2:end,:); x3(2:end,:)];

    clear Flux
    clear J_O2
        for j = 1:length(t)
        [Flux, Fluxbox(j,:)] = Mito_Flux(x(j,:),param1,boxVmax,delta_sf);
        J_O2(j) = Flux(6)/2; % J_O2 in mol / sec / L mito         %% *Rm_cell*Rc_t/(1000*rho)*60*1e6; % (umol O2) min-1 (g tissue)-1 . %% Flux(4)/2*1.8957e5
  
        end
        clear J_electrode

% Accounting for electrode response time
    [t_O2,J_electrode] = ode15s(@dxdt_electrode, t, 0, [], t, J_O2);  
    tss=find(t>90,1);
    
    J_O2_el_n(ii) = J_electrode(tss); % max(J_electrode)
     
    end
    
    sim_JO2_nocarn(kk,:)=J_O2_el_n;
end
toc

save('Variability_Run')