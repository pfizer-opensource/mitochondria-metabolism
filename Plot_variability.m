load('Variability_Run')
%SA_isolated_mito_experiment % Raw code to produce variability run

%% Plot without L-Carnitine
JO2_nocarn=mean(sim_JO2_nocarn,1);
JO2_nocarnstd=std(sim_JO2_nocarn,1);
 Obs_JO2n = [ 9.8278E-05; 0.00012433; 0.000294761; 0.000652426; 0.000629419; 0.000659878; 0.000591224];
 obserrn=[1.22185E-05;1.06737E-05;5.96623E-05;0.000112975;0.000133269; 0.000170832; 0.000130777];
    Fit_Jn = flip(JO2_nocarn);
    Jabsn=[Obs_JO2n,Fit_Jn'];
    xbar=[1,2,3,4,5,6,7];
    figure;
    hl=bar(xbar,Jabsn,'grouped');
    hold on
    labels={'C4','C6','C8','C10','C12','C14','C16'};
    set(gca,'XTick',1:7,'XTickLabel',labels)
    set(gca,'Fontsize',20)
    
    set(hl(1), 'FaceColor', 'k')
    set(hl(2), 'FaceColor', [0.5 0.5 0.5])
    ylabel('VO_2 max (mol/s/L mito)')
    % title('Without L-Carnitine in the cytosol');
    ylim([0 1e-3])
  % Add errorbar
    ngroups = size(Jabsn, 1);
    nbars = size(Jabsn, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    
        for m= 1:2:nbars
        xpos=(1:ngroups) - (groupwidth/2) + (2*m-1) * (groupwidth / (2*nbars));
        hold on
        errorbar(xpos,Jabsn(:,m),obserrn,obserrn,'color',[0.5 0.5 0.5],'linestyle','none')
        end
        hold on
        for m= 2:1:nbars
        xpos=(1:ngroups) - (groupwidth/2) + (2*m-1) * (groupwidth / (2*nbars));
        hold on
        errorbar(xpos,Jabsn(:,m),JO2_nocarnstd,JO2_nocarnstd,'color',[0.5 0.5 0.5],'linestyle','none')
        end
        lgd1=legend('Observed JO_2','Model Simulated JO_2','Standard Deviation');
        lgd1.Location = 'northwest'; lgd1.FontSize=12;
        xlabel('C_n  Acyl-Carnitine');

%% Plot with L-Carnitine
        JO2_carn=mean(sim_JO2_Carn,1);
        JO2_carnstd=std(sim_JO2_Carn,1);

        Obs_JO2 =[0.000477126; 0.000554358;  0.000547048;  0.000601049;  0.000272728;  0.000123818; 8.07258E-05;];
        Obs_JO2 =flip(Obs_JO2); 
        obserr=[5.11352E-05; 1.68106E-05; 2.39637E-05; 4.97114E-05; 0.000112614; 0.000157951; 9.70213E-05];

        Fit_J=flip(JO2_carn);
        Jabsn=[Obs_JO2,Fit_J'];
        
        xbar=[1,2,3,4,5,6,7];
        figure;
        h=bar(xbar,Jabsn,'grouped');
        hold on
        labels={'C4','C6','C8','C10','C12','C14','C16'};
        set(gca,'XTick',1:7,'XTickLabel',labels)
        set(gca,'Fontsize',20)
        set(h(1), 'FaceColor', 'k')
        set(h(2), 'FaceColor', [0.5 0.5 0.5])
        ylabel('VO_2 max (mol/s/L mito)')
        %title('With L-Carnitine in the cytosol');
        ylim([0 1e-3])
    
    % Add errorbar
        for m= 1:2:nbars
        xpos=(1:ngroups) - (groupwidth/2) + (2*m-1) * (groupwidth / (2*nbars));
        hold on
        errorbar(xpos,Jabsn(:,m),obserr,obserr,'color',[0.5 0.5 0.5],'linestyle','none')
        end
        hold on
        for m= 2:1:nbars
        xpos=(1:ngroups) - (groupwidth/2) + (2*m-1) * (groupwidth / (2*nbars));
        hold on
        errorbar(xpos,Jabsn(:,m),JO2_carnstd,JO2_carnstd,'color',[0.5 0.5 0.5],'linestyle','none')
        end
     	
     	lgd= legend('Observed JO_2','Model Simulated JO_2','Standard Deviation')
    	lgd.Location = 'northwest'; lgd.FontSize=12;
    	xlabel('C_n  Acyl-Carnitine');
    