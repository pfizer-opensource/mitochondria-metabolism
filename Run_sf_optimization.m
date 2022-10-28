%function [Oss_out,Fit_JO2]=Run_optimization

delta_sf =[0.006; 0.003; 1; 0.5; 0.5;]; 

%boxVmax = [7e-9; 7e-9; 6.5167e-9; 1.3333e-10; 1.6667e-10;  1.35e-9;  1.35e-9;  6.0e-8;  1.6667e-8;  6.2833e-9;  4.7333e-8; 2e-10;]; % Bakker Vmax values (M/s)

%% Define bounds

p_lower=delta_sf.*0.5;
p_upper=delta_sf.*2;


%% fmincon
options = optimoptions('fmincon','Display','iter');
options.Algorithm = 'interior-point';
[sf_out,Fval,~,~,~,~,Hess]=fmincon(@(delta_sf)Objfun(delta_sf),delta_sf,[],[],[],[],p_lower,p_upper,[],options);

% Save Results
save('sf_optimization')
save('sf_out','sf_out')



