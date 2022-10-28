function [J_PDH_matrix]=PDHphosdephos(acetylcoA_matrix,ATP_matrix,ADP_matrix,NADH_matrix,COAS_matrix,NAD_matrix,pyruvate_matrix,CO2tot_matrix,Vmf_PDH,Keq_PDH_matrix,ai2,ai1,ref_pyr,Varstruc)
% ACCOA_x,ATP_x,ADP_x,NADH_x,COASH_x,NAD_x,PYR_x,CO2tot,Vmf,Keq_pdh,ai1,ai2
% PDHphosdephos(ACCOA_x,ATP_x,ADP_x,NADH_x,COASH_x)
% This function describes a simple way to incorporate activation and
% inactivation by phosphorylation of PDH (DAB)

%% Derivation from inhibition kinetics
% Where J is flux, V is Vmax, (Vr is reverse), Ka is association constant, 
% Ki is the inhibition constant, and A is concentration of substrate
% J+ = V+ / 1+ ka/(1-a)
% J- = V- /1 + ki/a
% J+ = J-
% V-/V+ = (1 + Ki/a)/(1+ka/(1-a))
% **Simplification:  ki = ka = 1/2 
% Vr = V- / V+
% aa = (3Vr-1 - sqrt(9 Vr^2 - 14Vr +9))/4Vr-4
% aa= 1/2 for Vr = 1
%


%      ref conc. 1.944e-3 was determined by starting with ic and manually titrating up
rat=(acetylcoA_matrix/COAS_matrix)*(NADH_matrix/NAD_matrix)*(ATP_matrix/ADP_matrix)*(ref_pyr/(pyruvate_matrix+1e-20));
if(rat>0.9999)&&(rat<1.0001)
  aa=0.5;
elseif Varstruc.PDH_UR == 1
 aa = 1;
else
  aa=Varstruc.PDH_DR*(3*rat-1-sqrt(9*(rat^2)-14*rat+9))/(4*rat-4); 
end

save('aa','aa')

a=pyruvate_matrix;
b=COAS_matrix;
c=NAD_matrix;
p=CO2tot_matrix;
q=acetylcoA_matrix;
r=NADH_matrix ;                                                                         
KmA=38.3e-6;
KmB=9.9e-6;
KmC=60.7e-6;
KiACCOA=40.2e-6;
KiNADH=40.0e-6;



J_PDH_matrix=aa*Vmf_PDH*(a*b*c-p*q*r/Keq_PDH_matrix)/(KmC*ai2*a*b+KmB*ai1*a*c+KmA*b*c+a*b*c);

