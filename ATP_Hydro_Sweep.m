function [x_st4, Flux_st4,dGATP,PCRATP,pdhcs,pi_c,MVO2,ADratio]=ATP_Hydro_Sweep(xo,param,Varstruc)

%% Globals
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
       iC14Carn_cy iC14Carn_m iC14CoA_m iC14EnoylCoA_m iC14OHCoA_34e4r5m...    
       iC14KetoCoA_m iC12Carn_cy iC12Carn_m iC12CoA_m iC12EnoylCoA_m iC12OHCoA_m...   
       iC12KetoCoA_m iC10Carn_cy iC10Carn_m iC10CoA_m iC10EnoylCoA_m iC10OHCoA_m...  
       iC10KetoCoA_m iC8Carn_cy iC8Carn_m iC8CoA_m iC8EnoylCoA_m iC8OHCoA_m...  
       iC8KetoCoA_m iC6Carn_cy iC6Carn_m iC6CoA_m iC6EnoylCoA_m iC6OHCoA_m...
       iC6KetoCoA_m iC4Carn_cy iC4Carn_m iC4CoA_m...
       iC4EnoylCoA_m iC4OHCoA_m iC4KetoCoA_m iAcetylCoAMAT iFADH_m iNADHm iCoAMAT...
       iC16AcylCoACYT iC16Carn_cy;
%%
colors=['k','b','r','g','m'];
Vmito=2.8820e-1*Varstruc.mito_rho; % Mitochondria Volume modulatd by rho to simulate mito dysfunction
Vcyto= 0.94-Vmito;  			   % Cytoplasmic volume (less than 1) 
W_x = 6.5138e-1;				   %
W_i = 7.2376e-2;
W_c = 8.4251e-1;
Rm_cell = Vmito;
Rc_cell = Vcyto; 
Rc_t =  7.3078e-01; 
rho = 1.0500;


Np=linspace(0.36e-3,1.3e-3,10);  % ATP hydrolysis rates (M/s)
NP=Np*1000;						 % Convert hydrolysis rate to mM/s
param1=param; 					 % Define which parameters to use
param1(44)=Varstruc.NADtot_p44;	 %
iflag=100;						 % in vivo flag
options=odeset('RelTol',1e-8,'AbsTol',1e-12);  % Set toleranc
xo1=xo;							 % Define initial Conditions

tic
for j = 1: numel(Np)
   
   x_AtC(j) = Np(j)  ;%(j-1)*(0.8e-3)/10; % mol s^{-1} (l cyto)^{-1}, ATP consumption rate
   clear xnew
   clear t4
   clear Flux4
   [t4,xnew] = ode15s(@Cell_dxdT,[0 2000],xo1,options,param1,x_AtC(j),iflag,Varstruc);
   x_st4(j,:)= xnew(end,:); 
   Flux_st4(j,:) = Cell_Flux(xnew(end,:),param1,x_AtC(j),Varstruc);
   ATPobs(j) = W_c*Rc_cell*xnew(end,iATP_c) + W_x*Rm_cell*xnew(end,iATP_x) + W_i*Rm_cell*xnew(end,iATP_i);
   PCRobs(j) = W_c*Rc_cell*xnew(end,iPCr_c);
   ATP_c(j) = xnew(end,iATP_c);
   ADP_c(j) = xnew(end,iADP_c);
   PI_c(j) = xnew(end,iPI_c);
   ADPobs(j)=W_c*Rc_cell*xnew(end,iADP_c)+W_x*Rm_cell*xnew(end,iADP_x)+ W_i*Rm_cell*xnew(end,iADP_i);
   AMPobs(j)=W_c*Rc_cell*xnew(end,iAMP_c)+W_x*Rm_cell*xnew(end,iAMP_x)+ W_i*Rm_cell*xnew(end,iAMP_i);
   MVO2(j) = (Flux_st4(j,6)./2).*Rm_cell*Rc_t/(1000*rho)*60*1e6; % Convert units from M/s to (umol O2) min-1 (g tissue)-1

end
 
pi_c = x_st4(:,iPI_c); 
dGATP = Flux_st4(:,46);
PCRATP = PCRobs./ATPobs;

ADratio.pcratp=PCRATP;
ADratio.atpadp=ATPobs./ADPobs;
ADratio.atpamp=ATPobs./AMPobs;
ADratio.amp=AMPobs;
ADratio.atp=ATPobs;
ADratio.adp=ADPobs;

pdh = Flux_st4(:,12);
cs = Flux_st4(:,13);
pdhcs=pdh./cs;

