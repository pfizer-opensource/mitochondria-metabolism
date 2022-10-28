function [Oss_out]=Objfun(delta_sf)
%% Define Globals
boxVmax = [7e-9; 7e-9; 6.5167e-9; 1.3333e-10; 1.6667e-10;  1.35e-9;  1.35e-9;  6.0e-8;  1.6667e-8;  6.2833e-9;  4.7333e-8; 2e-10;]; % Bakker Vmax values (M/s)

setup;
param(39) = 250; % leak parameter
param(1)=param(1)/100; % Calibratd for pdh model substrate selection    PDH . - 2.05e-5 M/s
param(2)= param(2)/100;% Calibratd for pdh model substrate selction   CS . - 9.82e-4 M/s
param(44)= 5.6e-5; % ki of CnAcyl carnitines in the cytosol for CACT 56 uM
param(45)= 2.0e-4; % ki of CACT for carnitine in cytosol


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
SF=10^6; % M
xo(105) = 0;
xo(106) = 0.46/SF;  % FADH2 Matrix 
xo(108) =0;         % CoASH Rule 
xo(110)=25/SF;      % 25 uM Carnitine in cytosol
xo(111)=25/SF;      % 25 uM Carnitine in the intermembrane space

%% Run IsoMito experiment with L-Carnitine
 
xo(iACCOA_x) = 70e-6; 
xo(iNADH_x) = 16/SF; 
xo([42:62 iATP_c iADP_c iAMP_c iATP_i iADP_i iAMP_i]) = 0;
[t0,x0] = ode15s(@Mito_dXdT,[0 60],xo,options,param,boxVmax,delta_sf);
 
 
% Do state-1 experiment with experimental buffer levels % 1.6 seconds 
 
xo1 = x0(end,:);
xo1(iPI_c) = 5.0e-3;
xo1(iH_c) = 10.^(-7.2);
[t1,x1] = ode15s(@Mito_dXdT,[0 60],xo1,options,param,boxVmax,delta_sf);
 
Carnind = [63,69,75,81,87,93,99];
for i = 1:numel(Carnind)

% Add substrates and run state-2 simulation 
xo2 = x1(end,:);
xo2(iMAL_c) = 1.0e-3;
xo2(Carnind(i)) = 20e-6;
[t2,x2] = ode15s(@Mito_dXdT,[0 60],xo2,options,param,boxVmax,delta_sf);
 
 
% Add ADP for state 3 simulation 
xo3 = x2(end,:);
xo3(iADP_c) = 0.375e-3;
[t3,x3] = ode15s(@Mito_dXdT,[61 120],xo3,options,param,boxVmax,delta_sf);
 
% Combine simulation timecourses
clear t
clear x
t = [t2(2:end); t3(2:end)]; 
x = [x2(2:end,:); x3(2:end,:)];
 
clear J_O2
for j = 1:length(t)
  Flux = Mito_Flux(x(j,:),param,boxVmax,delta_sf);
  J_O2(j) = Flux(6)/2; % J_O2 in mol / sec / L mito
end
% Accounting for electrode response time
 [t_O2,J_electrode] = ode15s(@dxdt_electrode, t, 0, [], t, J_O2);  
  tss=find(t>=90,1);
 J_O2_el(i) =  J_electrode(tss);
 
end




%% Run IsoMito experiment without L-Carnitine

xo(110) = 0; 
xo(iACCOA_x) = 70e-6; 
xo(iNADH_x) = 16/SF; 
xo([42:62 iATP_c iADP_c iAMP_c iATP_i iADP_i iAMP_i]) = 0;
[t0n,x0n] = ode15s(@Mito_dXdT,[0 60],xo,options,param,boxVmax,delta_sf);


% Do state-1 experiment with experimental buffer levels  

xo1n = x0n(end,:);
xo1n(iPI_c) = 5.0e-3;
xo1n(iH_c) = 10.^(-7.2);
[t1n,x1n] = ode15s(@Mito_dXdT,[0 60],xo1n,options,param,boxVmax,delta_sf);

Carnind = [63,69,75,81,87,93,99];
for k = 1:numel(Carnind)
    k
% Add substrates and run state-2 simulation 
xo2n = x1n(end,:);
xo2n(iMAL_c) = 1.0e-3;
xo2n(Carnind(k)) = 20e-6;
[t2n,x2n] = ode15s(@Mito_dXdT,[0 60],xo2n,options,param,boxVmax,delta_sf);


% Add ADP for state 3 simulation 

xo3n = x2n(end,:);
xo3n(iADP_c) = 0.375e-3;
[t3n,x3n] = ode15s(@Mito_dXdT,[61 120],xo3n,options,param,boxVmax,delta_sf);

% combine simulation timecourses
clear tn
clear xn
tn = [t2n(2:end); t3n(2:end)]; 
xn = [x2n(2:end,:); x3n(2:end,:)];

clear Fluxn
clear J_O2n
for m = 1:length(tn)
  Fluxn = Mito_Flux(xn(m,:),param,boxVmax,delta_sf);
  J_O2n(m) = Fluxn(6)/2; % J_O2 in mol / sec / L mito
end

clear J_electroden
% Accounting for electrode response time
 [t_O2n,J_electroden] = ode15s(@dxdt_electrode, tn, 0, [], tn, J_O2n);  
 tss=find(tn>=90,1);
 
 J_O2_eln(k) =  J_electroden(tss);
    
end


Fit_JO2=[J_O2_el'; J_O2_eln';];

            %C16Carn; C14Carn;       C12 Carn;     C10Carn;       C8Carn; %%C6 Carn;       C4Carn;          Without CArn 
Obs_JO2=[0.000477126; 0.000554358;  0.000547048;  0.000601049;  0.000272728;  0.000123818; 8.07258E-05; 0.000591224; 0.000659878; 0.000629419; 0.000652426;0.000294761; 0.00012433; 9.8278E-05;]; 

Oss=(Fit_JO2-Obs_JO2).^2 ;
 W =Oss.*(1./(Fit_JO2).^2); % Weighting by 1/Y_pred^2
 Oss_out=sum(W);
 


end
