function Run_substrate_select(xo,param,Varstruc,x_AtC)

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
%colors=['k','b','r','g','m'];
W_x = 6.5138e-1;
W_i = 7.2376e-2;
W_c = 8.4251e-1;
Rm_cell = 2.8820e-1;
Rc_cell = 6.8010e-1; 
Rc_t =  7.3078e-01; 
rho = 1.0500;

param1=param; 
param1(44)=Varstruc.NADtot_p44;
iflag=100;
options=odeset('RelTol',1e-8,'AbsTol',1e-12);
xo1=xo;

tic
Parmvary=[0.1,1,10,100,1000];
         
        for i = 1:numel(Parmvary)
            
             xo2=xo1;
             xo2(iPYR_c)=(xo2(iPYR_c)*Parmvary(i)); %xo2(iPYR_c)+
             for m = 1: numel(Parmvary)
                 
                 xo2(109)=xo1(109);
                 xo2(109)=(xo2(109)*Parmvary(m)); %xo2(109)+

                clear xnew2
                clear t3
                [t3,xnew2] = ode15s(@Cell_dxdT,[0 2000],xo2,options,param1,x_AtC,iflag,Varstruc);

                    for k = 1:length(t3)
                    Flux2(k,:) = Cell_Flux(xnew2(k,:),param1,x_AtC,Varstruc);
                    end


                    if i==1
                        f=m;
                    else
                        f=(m+(5*(i-1)));
                    end
                    x_st3(f,:)= xnew2(end,:);
                    Flux_st3(f,:) = Cell_Flux(xnew2(end,:),param1,x_AtC,Varstruc);

             end 
        end
        
        
        % Plots
        PDH_J=Flux_st3(:,12);
        CS_J=Flux_st3(:,13);
        PDH_CS=PDH_J./CS_J;

        PCoA=x_st3(:,109);
        Pyr=x_st3(:,iPYR_c);

        Z=reshape(PDH_CS,[5,5]);
        Pyr_s=Pyr([1:5:25]);
        PCoA_s=PCoA([1:5]);

        figure(4);
        surf(Pyr_s,PCoA_s,Z)
        xlabel('Pyruvate Conc. (M)')
        ylabel('Palmitoyl-CoA Conc. (M)')
        zlabel('PDH/CS Ratio')
        set(gca,'Xscale','log','Yscale','log')
        zlim([0 1.2])
        zticks([0.2,0.4,0.6,0.8,1,1.2]);
        yticks([0.0001, 0.01,  1]);
        xticks([0.00001, 0.001]);
        if x_AtC  == 0.5e-3
        title({'Substrate Selection at Rest';'ATP Hydrolysis Rate (0.5 mM/s)';})
        elseif x_AtC  == 1.0e-3
        title({'Substrate Selection Under Moderate Exercise';'ATP Hydrolysis Rate (1.0 mM/s)';})
        end
        
        
        pcr = W_c*Rc_cell*x_st3(:,iPCr_c);
        atp = W_c*Rc_cell*x_st3(:,iATP_c) + W_x*Rm_cell*x_st3(:,iATP_x) + W_i*Rm_cell*x_st3(:,iATP_i);

        pcratp = pcr./atp;
        zz=reshape(pcratp,[5,5]);
        

end
