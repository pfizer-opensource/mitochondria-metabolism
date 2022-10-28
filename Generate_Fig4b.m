function [sim_pop]=Generate_Fig4b(Ind,X_2,xo1,param1,Varstruc)
% Simulate Energetics and Substrate Selection from rest to exercise
% Ind is the Index for each paramter sampled, X_2 is the latin hypercube
% of parameter sets based on sig variability, xo1 is the initial conditions
nInd=numel(Ind);
for k=1:1:length(X_2(:,1))
        
        for m=1:nInd
        param1(Ind(m)) =X_2(k,m); 
        end
        
        for n=1 
        Varstruc.boxVmax(12)=X_2(k,n+nInd); 
        end
        
        for nn=1
        Varstruc.Kmcpt1CarCYT =X_2(k,nn+nInd+1);
        end
 
     [~,~,dGATP,PCRATP,pdhcs,pi,MVO2]=ATP_Hydro_Sweep(xo1,param1,Varstruc);
     sim_pop.pcr_atph(1,:)= dGATP;
     sim_pop.pcr_atph(k+1,:)=PCRATP;
     sim_pop.pdhcsh(1,:)=dGATP;
     sim_pop.pdhcsh(k+1,:)=pdhcs;
     sim_pop.pih(1,:)=dGATP;
     sim_pop.pih(k+1,:)=pi*1000; % convert to mM 
     sim_pop.MVO2(1,:)= dGATP;
     sim_pop.MVO2(k+1,:)=MVO2;
end
