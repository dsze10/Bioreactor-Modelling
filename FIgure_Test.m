close all
%Figures using Plateau Function [y] = BioReactor_model(k,yo,ym,t)

%Theoretical Maximum DNA (Minf)
Minf = figure('visible','on');

Tinf = 0:100;
TotalDNA = 100;
ExtractableDNA = 60;
Y0 =0;
K = 0.05;
    
Texp = 0:30;
Yexp = BioReactor_model(K,Y0,ExtractableDNA,Texp);

plot(Tinf,BioReactor_model(K,Y0,TotalDNA,Tinf),...
    Tinf,BioReactor_model(K,Y0,ExtractableDNA,Tinf),...
    Texp,Yexp,...
        'LineWidth',3,'LineStyle','-');
 yline(TotalDNA,'--',{'Total DNA in Tissue'});
 yline(ExtractableDNA,'--',{'DNA extractable by Reagent'});
 yline( Yexp(end),'--',{'DNA extracted in Exp'});
 
 set(gca,'FontSize',16);
    ylabel('DNA','FontSize',20);
    xlabel('Time','FontSize',20);
    title(strcat('Theoretical Max DNA (M\infty)'),'FontSize',24);
    ylim([0 TotalDNA+10]);
    xlim([0 inf]);
    
saveas(Minf,char(strcat('Supp/','Minf.pdf')));

% % 
%% Param estimations Fig 7


%Theoretical Maximum DNA (Minf)
Fig7PE = figure('visible','on');
Params = array2table(zeros(4,length(Params)), 'VariableNames',Params);
Params{1,:} = [0.971435676,0.01071063,0.609174714,0.061296883];
Params{2,:} = [0.896175012,0.012157603,0.636214571,0.059381464];
Params{3,:} = [0.92598316,0.011436507,0.639378824,0.057026406];
Params{4,:} = [0.94409352,0.011166094,0.827617833,0.031083404];



Texp = EstData.Time;
Yexp = EstData.Data;
TotalDNA; 

plot(Texp,Yexp,...
        'LineWidth',3,'LineStyle','-');
    hold on
    
    for c= 1:size(Params,1)
        
       Y_n=BioReactor_model(,y0n,ymn,0:tndur);
       plot(tn0:tnend,Yexp,...
        'LineWidth',3,'LineStyle','-'); 
       Y_d=BioReactor_model(kn,y0n,ymn,0:tddur);
       plot(td0:tdend,Yexp,...
        'LineWidth',3,'LineStyle','-'); 
    end
    

 
 set(gca,'FontSize',16);
    ylabel('DNA','FontSize',20);
    xlabel('Time','FontSize',20);
    title(strcat('Parameter Estimation'),'FontSize',24);
    ylim([0 TotalDNA+10]);
    xlim([0 inf]);
    
saveas(Minf,char(strcat('Supp/','Minf.pdf')));
%% 

%Ideal Optimization of Model 

Opti = figure('visible','on');
TotalDNA = 100;
Y0 =0;
MaxExtract_N = 0.8;

% Equal Time High Kd
Tn = 0:60;
Kn = 0.01;
YMn = TotalDNA*MaxExtract_N;

Yn = BioReactor_model(Kn,Y0,YMn,Tn);

Td = 0:60;
Kd = 0.05; 
Ymd = TotalDNA-Yn(end);

Yd = BioReactor_model(Kd,Y0,Ymd,Td);

plot(Tn,Yn,...
        Td+Tn(end),Yd+Yn(end),...
        'LineWidth',3,'LineStyle','-','Color','b');
    
% Optimized time using inverse
Opti = 0.8

Tn = 0:60;
Kn = 0.01;
YMn = TotalDNA*MaxExtract_N;

Yn = BioReactor_model(Kn,Y0,YMn,Tn);

Td = 0:60;
Kd = 0.05; 
Ymd = TotalDNA-Yn(end);

Yd = BioReactor_model(Kd,Y0,Ymd,Td);

plot(Tn,Yn,...
        Td+Tn(end),Yd+Yn(end),...
        'LineWidth',3,'LineStyle','-','Color','b');
    
set(gca,'FontSize',16);
    ylabel('DNA','FontSize',20);
    xlabel('Time','FontSize',20);
    title(strcat('Total DNA extracted'),'FontSize',24);
    ylim([0 TotalDNA+10]);
    xlim([0 inf]);

    %% 
    createfigure()
    
 figs= {'IT1_Fig7PE_0.1_0.01.fig','IT2_Fig7PE_0.1_0.01.fig','IT3_Fig7PE_0.1_0.01.fig','IT4_Fig7PE_0.1_0.01.fig'};
x={[];[];[];[]};
y={[];[];[];[]};
z={[];[];[];[]};

 for fig = figs
    
     openfig(fig{1}); 

    a = get(gca,'Children');
    xdata = get(a, 'XData');
    ydata = get(a, 'YData');
    zdata = get(a, 'ZData');
    
  x{:,1} = ;
  test = [x{:,1},xdata{:,1}];

 end

 %% 
 pre = openfig('Pre_Fig3sim.fig');
 post= openfig('POST_Fig3sim.fig');
 
L = findobj(pre,'type','line');
copyobj(L(2),findobj(post,'type','axes'));
%fig7
 pre = openfig('Pre_Fig7sim.fig');
 post= openfig('POST_Fig7sim.fig');
 
L = findobj(pre,'type','line');
copyobj(L(2),findobj(post,'type','axes'));
%%
close all
%Figures using Plateau Function [y] = BioReactor_model(k,yo,ym,t)

%Theoretical Maximum DNA (Minf)
Minf = figure('visible','on');

Tinf = 0:100;
TotalDNA = 100;
ExtractableDNA = 60;
Y0 =0;
K = 0.05;
    
Texp = 0:30;
Yexp = BioReactor_model(K,Y0,ExtractableDNA,Texp);

plot(Tinf,BioReactor_model(K,Y0,TotalDNA,Tinf),...
        'LineWidth',3,'LineStyle','-');
 yline(TotalDNA,'--',{'YM'});
 yline(0,'--',{'Y0'});
 
 set(gca,'FontSize',16);
    ylabel('','FontSize',20);
    xlabel('Time','FontSize',20);
    title(strcat('Exponential Plateau Equation'),'FontSize',24);
    ylim([0 TotalDNA+10]);
    xlim([0 inf]);
    
saveas(Minf,char(strcat('Supp/','Minf.pdf')));

%%
%Theoretical Maximum DNA (Minf)
Minf = figure('visible','on');

Tinf = 0:100;
TotalDNA = 100;
ExtractableDNA = 60;
Y0 =0;
K = 0.05;
    
Texp = 0:30;
Yexp = BioReactor_model(K,Y0,ExtractableDNA,Texp);

plot(Tinf,BioReactor_model(K,Y0,TotalDNA,Tinf),...
    Tinf,BioReactor_model(K,Y0,ExtractableDNA,Tinf),...
    'LineWidth',2)
 yline(TotalDNA,'--',{'Total DNA in Tissue'},'FontSize',16);
 yline(ExtractableDNA,'--',{'DNA extractable by Reagent'},'FontSize',16);
 
 set(gca,'FontSize',16);
    ylabel('DNA','FontSize',20);
    xlabel('Time','FontSize',20);
    title(strcat('Theoretical Maximum DNA (YM)'),'FontSize',24);
    ylim([0 TotalDNA+10]);
    xlim([0 inf]);


