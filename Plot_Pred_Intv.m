function Plot_Pred_Intv(BW_tpl,color_flag,fig_num)

colors=['k','b','r','g','m'];
yla={'PCr/ATP Ratio','J_P_D_H/J_C_S','Inorganic Phosphate (Pi, mM)','PCr/ATP Ratio','J_P_D_H/J_C_S Ratio','Pi (mM)'};
shade(1).struc(:)=[0.5 0.5 0.5];
shade(2).struc(:)=[0.5843 0.8157 0.9882]; %light blue  %cyan= [0.2 0.8 0.8];
shade(3).struc(:)=[1 0.2 0.2];
shade(4).struc(:)=[0 0.6 0.3];
shade(5).struc(:)=[1 0 1];
% ATP Hydrolysis rate
if fig_num == 1 || fig_num == 2 || fig_num == 3
NP = BW_tpl(1,:);
else
Np=linspace(0.36e-3,1.3e-3,5); 
NP=Np*1000;
end


% Plotting Mean and 95% Prediction Interval
  BWoutmean=mean(BW_tpl([2:end],:));
  BWoutstd=std(BW_tpl([2:end],:));
  N = numel(BW_tpl([2:end],1)); 
  BWoutSEM=BWoutstd/sqrt(N); %Can be used to calculate 95% CI of the mean
  CI95 = tinv([0.05 0.95], N-1); 
  yCI95 = bsxfun(@times, BWoutstd, CI95(:));   % Calculate 95 CI Prediction interval
  y95_lower=yCI95(1,:)+BWoutmean;
  y95_upper=yCI95(2,:)+BWoutmean;
  x=NP; %BW_tpl(1,:);
  x2 = [x, fliplr(x)];
  inBetween = [y95_lower, fliplr(y95_upper)];
  
  figure(fig_num);
  hold on
  plot(NP,BWoutmean,colors(color_flag),'linewidth',2) %BW_tpl(1,:)
  fill(x2, inBetween, shade(color_flag).struc,'linestyle','none');
  alpha(0.25) %Make fill transparent
  ylabel(yla(fig_num)); 
  if fig_num == 1 || fig_num == 2 || fig_num == 3
  xlabel('\DeltaG_A_T_P');
  else
      xlabel('ATP Hydrolysis (mM/s)');
  end
      
  axis auto;

  