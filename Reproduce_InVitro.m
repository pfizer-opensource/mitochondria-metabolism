%% Reproducing in vitro data from:
% K. H. Fisher-Wellman, M. T. Davidson, T. M. Narowski, C.-T. Lin, T. R. Koves and D. M. Muoio
% Cell reports 2018 Vol. 24 Issue 13 Pages 3593-3606.e10
% Accession Number: 30257218 DOI: 10.1016/j.celrep.2018.08.091

% Publication Figure 3A
StateType = 3 ;   % Initialize the mitochondrial oxidative state to state 3
ExpType   = 7 ;   % Specify experimental initial conditions
fig_flag  = 1 ;   % Specify which plots to generate

InVitroExpCond(fig_flag,StateType,ExpType) % Run function

% Publication Figure 3B
fig_flag =2 ;
InVitroExpCond(fig_flag,StateType,ExpType)

