% Modeling and Programming Methods to Streamline Biomaterial Development
% Daphne Sze 
% November 2022

% Matlab script for plotting model with YM
%% Run with Model with Fitted Ym and Total DNA YM
clc
close all

% Fixed Constants
t0=0;
ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
Reactor_Volume = 25; %ml
hD2W =0.3; %Dry to wet ratio
lD2W =0.2; 
highDNA = 14.5; %DNA concentration ug/mg dry weight
lowDNA = 9;

%Initialize array for figuers 
pobj_ym = gobjects(2,size(Fitdat,2)-1);

%For all datasets
for c= Fitdat(1,2:end)
    %Get column for dataset
    col = find(strcmp(Fitdat(1,:),c));
    
    %Get data from table
    NaDeoDat = Fitdat{5,col};
    NaDeoParam = Fitdat{6,col};
    DNaseDat = Fitdat{7,col};
    DNaseParam = Fitdat{8,col};
   
    %Get Mass and calculate total DNA ranges
    Mass_ECM_wet = Mass{strcmp(Mass.Fig, c), 3}; %mg
    lMass_ECM_dry = Mass_ECM_wet*lD2W; 
    hMass_ECM_dry = Mass_ECM_wet*hD2W; 
    lltotal_DNA =  lowDNA*lMass_ECM_dry;
    hltotal_DNA =  highDNA*lMass_ECM_dry;
    lhtotal_DNA =  lowDNA*hMass_ECM_dry;
    hhtotal_DNA =  highDNA*hMass_ECM_dry; 

 %Get plotting Data
    %FitDat row 6:NaDeo ; row 8: Dnase
    n_duration= NaDeoParam{strcmp(NaDeoParam.Param, 'Dur'), 2}; %Duration of the Sodium Doexycholate Step in min 
    d_duration=DNaseParam{strcmp(DNaseParam.Param, 'Dur'), 2}; %%Duration of the DNAse Step in min
       
%For NaDeoxy
    tnstop=t0+n_duration; %min
    kn= NaDeoParam{strcmp(NaDeoParam.Param, 'k'), 2}; %Overall rate constant min^-1
    y0n=0; % Initial DNA concentration [ug/ml] after PBS wash
    
    %Time
    T_n=t0:n_duration;
    
    %Calculate YM and model results
    % 9ug/mg, 20% D:W
    llymn=(lltotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    llY_n=BioReactor_model(kn,y0n,llymn,T_n)*ConversionFactor* Reactor_Volume;
    % 14.5 ug/mg, 20% D:W
    hlymn=(hltotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    hlY_n=BioReactor_model(kn,y0n,hlymn,T_n)*ConversionFactor* Reactor_Volume;
    % 9ug/mg, 30% D:W
    lhymn=(lhtotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    lhY_n=BioReactor_model(kn,y0n,lhymn,T_n)*ConversionFactor* Reactor_Volume;
    % 14.5 ug/mg, 30% D:W
    hhymn=(hhtotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    hhY_n=BioReactor_model(kn,y0n,hhymn,T_n)*ConversionFactor* Reactor_Volume;
    %Prism Fit
    ymn=NaDeoParam{strcmp(NaDeoParam.Param, 'YM'), 2};
    Y_n=BioReactor_model(kn,y0n,ymn,T_n)*ConversionFactor* Reactor_Volume;
    
%Get Twash
    twash = DNaseParam{strcmp(DNaseParam.Param, 'twash'), 2};

%For DNAse
    tdstop=tnstop+twash+d_duration; %min
    kd= DNaseParam{strcmp(DNaseParam.Param, 'k'), 2}; %Overall rate constant min^-1
    y0d=0; % Initial DNA concentration [ug/ml] after PBS wash
    
    %Time 
    T_d=t0:d_duration;
    
    %Calculate YM and model results
    % 9ug/mg, 20% D:W
    llDNA_left=lltotal_DNA-llY_n(end);
    llymd=(llDNA_left/Reactor_Volume)/ConversionFactor;
    llY_d=BioReactor_model(kd,y0d,llymd,T_d)*ConversionFactor*Reactor_Volume;
    % 14.5 ug/mg, 20% D:W
    hlDNA_left=hltotal_DNA-hlY_n(end);
    hlymd=(hlDNA_left/Reactor_Volume)/ConversionFactor;
    hlY_d=BioReactor_model(kd,y0d,hlymd,T_d)*ConversionFactor*Reactor_Volume;
    % 9ug/mg, 30% D:W
    lhDNA_left=lhtotal_DNA-lhY_n(end);
    lhymd=(lhDNA_left/Reactor_Volume)/ConversionFactor;
    lhY_d=BioReactor_model(kd,y0d,lhymd,T_d)*ConversionFactor*Reactor_Volume;
    % 14.5 ug/mg, 30% D:W
    hhDNA_left=hhtotal_DNA-hhY_n(end);
    hhymd=(hhDNA_left/Reactor_Volume)/ConversionFactor;
    hhY_d=BioReactor_model(kd,y0d,hhymd,T_d)*ConversionFactor*Reactor_Volume;
    % Prism Fit
    ymd=DNaseParam{strcmp(DNaseParam.Param, 'YM'), 2};
    Y_d=BioReactor_model(kd,y0d,ymd,T_d)*ConversionFactor*Reactor_Volume;
    
 % Get Experimental Data for DNA content profiles
    cend = size(NaDeoDat,2);
    DNase = DNaseDat;
    DNase.('Time (min)')= DNase.('Time (min)')+twash+tnstop;
    ScatDat = [NaDeoDat;DNase];
    ScatDat.Mean =  mean(ScatDat{:,3:cend},2,'omitnan','double');
    ScatDat.STD = std(ScatDat{:,3:cend},0,2,'omitnan');
    ScatDat.SEM = ScatDat.STD./sqrt(sum(~isnan(ScatDat{:,3:cend}),2));
   
    
%Plot DNA profiles 
    x = ScatDat.('Time (min)')';
    y = ScatDat.Mean'*ConversionFactor* Reactor_Volume;
    yneg = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
    ypos = ScatDat.SEM'*ConversionFactor* Reactor_Volume;

    
    pobj_ym(1,col-1) = figure('visible','on');
    %Set plot parameters
    set(gca,'FontSize',16);
    ylabel('DNA [ug]','FontSize',20);
    xlabel('Time (min)','FontSize',20);
    title(strcat(c,' (',num2str(Mass_ECM_wet),' mg)'),'FontSize',24);
    ylim([0 inf]);
    xlim([0 inf]);
    
    %define x axis for wash
    con_x = [T_n(end) T_n(end)+twash];
    
    %Plot 9ug/mg, 20% D:W
    hold on
    plot(T_n,llY_n,...
        T_d+tnstop+twash,llY_d,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);
    hold on
    con_lly = [llY_n(end) llY_d(1)];
    cpll= plot(con_x,con_lly,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);

    %Plot 14.5 ug/mg, 20% D:W
    hold on
    plot(T_n,hlY_n,...
        T_d+tnstop+twash,hlY_d,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    hold on
    con_hly = [hlY_n(end) hlY_d(1)];
    cphl= plot(con_x,con_hly,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    
    %Plot 9ug/mg, 30% D:W
    hold on
    plot(T_n,lhY_n,...
        T_d+tnstop+twash,lhY_d,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    hold on
    con_lhy = [lhY_n(end) lhY_d(1)];
    cplh= plot(con_x,con_lhy,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);

    %Plot 14.5 ug/mg, 30% D:W
    hold on
    plot(T_n,hhY_n,...
        T_d+tnstop+twash,hhY_d,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    hold on
    con_hhy = [hhY_n(end) hhY_d(1)];
    cphh= plot(con_x,con_hhy,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    
    %Plot Prism Fit
    hold on
    plot(T_n,Y_n,...
        T_d+tnstop+twash,Y_d,...
        'LineWidth',3,'LineStyle','--','Color','b');
    hold on
    con_y = [Y_n(end) Y_d(1)];
    cp= plot(con_x,con_y,...
        'LineWidth',3,'LineStyle','--','Color','b');

    hold on
    %Fill regions between profiles
    patch([T_n fliplr(T_n)], [llY_n fliplr(hlY_n)], [147 197 207]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([con_x fliplr(con_x)], [con_lly fliplr(con_hly)], [147 197 207]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([T_d+tnstop+twash fliplr(T_d+tnstop+twash)], [llY_d fliplr(hlY_d)], [147 197 207]./255,'EdgeColor','none','FaceAlpha',.3)
    
    patch([T_n fliplr(T_n)], [hlY_n fliplr(lhY_n)], [156.5  206.5  194]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([con_x fliplr(con_x)], [con_hly fliplr(con_lhy)], [156.5  206.5  194]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([T_d+tnstop+twash fliplr(T_d+tnstop+twash)], [hlY_d fliplr(lhY_d)], [156.5  206.5  194]./255,'EdgeColor','none','FaceAlpha',.3)
    
    patch([T_n fliplr(T_n)], [lhY_n fliplr(hhY_n)], [166 216 181]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([con_x fliplr(con_x)], [con_lhy fliplr(con_hhy)], [166 216 181]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([T_d+tnstop+twash fliplr(T_d+tnstop+twash)], [lhY_d fliplr(hhY_d)], [166 216 181]./255,'EdgeColor','none','FaceAlpha',.3)
    
   %Error Bars
   hold on
   md= errorbar(x,y,yneg,ypos,'sb','MarkerSize',10,...
   'MarkerEdgeColor','black','MarkerFaceColor','black');

%Define Legend   
lgd = legend([cphh cplh cphl cpll cp md],'Model 14.5 ug/mg, 30% D:W','Model 9 ug/mg, 30% D:W','Model 14.5 ug/mg, 20% D:W','Model 9 ug/mg, 20% D:W','Model Fitted Ym','Experimental Data','FontSize',16,'Location','Northeast','NumColumns',1);
hold off

%Save Figures
saveas(pobj_ym(1,col-1),char(strcat('Figures/Ym/',c,'_1.pdf')));
saveas(pobj_ym(1,col-1),char(strcat('Figures/Ym/',c,'_1.fig')));


% Get Experimental Data for Total DNA extracted plots 

    cend = size(NaDeoDat,2);
    NaDeoend = mean(NaDeoDat{end,3:cend},'omitnan','double');
    DNase = DNaseDat;
    DNase.('Time (min)')= DNase.('Time (min)')+twash+tnstop;
    DNase{:,3:cend} = DNase{:,3:cend}+NaDeoend;
    ScatDat = [NaDeoDat;DNase];
    ScatDat.Mean =  mean(ScatDat{:,3:cend},2,'omitnan','double');
    ScatDat.STD = std(ScatDat{:,3:cend},0,2,'omitnan');
    ScatDat.SEM = ScatDat.STD./sqrt(sum(~isnan(ScatDat{:,3:cend}),2));

%Plot Total DNA plots
    x = ScatDat.('Time (min)')';
    y = ScatDat.Mean'*ConversionFactor* Reactor_Volume;
    yneg = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
    ypos = ScatDat.SEM'*ConversionFactor* Reactor_Volume;

    con_x = [T_n(end) T_n(end)+twash];
    
    pobj_ym(2,col-1) = figure('visible','on');
    
    %Set Graphing parameters
    set(gca,'FontSize',16);
    ylabel('DNA [ug]','FontSize',20);
    xlabel('Time (min)','FontSize',20);
    title(strcat(c,' (',num2str(Mass_ECM_wet),' mg)'),'FontSize',24);
    ylim([0 inf]);
    xlim([0 inf]);
    
    %Plot 9ug/mg, 20% D:W
    hold on
    plot(T_n,llY_n,...
        T_d+tnstop+twash,llY_d+llY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);
    hold on
    con_lly = [llY_n(end) llY_n(end)];
    cpll= plot(con_x,con_lly,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);

    %Plot 14.5 ug/mg, 20% D:W
    hold on
    plot(T_n,hlY_n,...
        T_d+tnstop+twash,hlY_d+hlY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    hold on
    con_hly = [hlY_n(end) hlY_n(end)];
    cphl= plot(con_x,con_hly,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    
    %Plot 9ug/mg, 30% D:W
    hold on
    plot(T_n,lhY_n,...
        T_d+tnstop+twash,lhY_d+lhY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    hold on
    con_lhy = [lhY_n(end) lhY_n(end)];
    cplh= plot(con_x,con_lhy,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    
    %Plot 14.5 ug/mg, 30% D:W
    hold on
    plot(T_n,hhY_n,...
        T_d+tnstop+twash,hhY_d+hhY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    hold on
    con_hhy = [hhY_n(end) hhY_n(end)];
    cphh= plot(con_x,con_hhy,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    
    %Plot Prism Fit
    hold on
    plot(T_n,Y_n,...
        T_d+tnstop+twash,Y_d+Y_n(end),...
        'LineWidth',3,'LineStyle','--','Color','b');
    hold on
    con_y = [Y_n(end) Y_n(end)];
    cp= plot(con_x,con_y,...
        'LineWidth',3,'LineStyle','--','Color','b');

    hold on
    %Fill regions between profiles
    patch([T_n fliplr(T_n)], [llY_n fliplr(hlY_n)], [147 197 207]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([con_x fliplr(con_x)], [con_lly fliplr(con_hly)], [147 197 207]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([T_d+tnstop+twash fliplr(T_d+tnstop+twash)], [llY_d+llY_n(end) fliplr(hlY_d+hlY_n(end))], [147 197 207]./255,'EdgeColor','none','FaceAlpha',.3)
    
    patch([T_n fliplr(T_n)], [hlY_n fliplr(lhY_n)], [156.5  206.5  194]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([con_x fliplr(con_x)], [con_hly fliplr(con_lhy)], [156.5  206.5  194]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([T_d+tnstop+twash fliplr(T_d+tnstop+twash)], [hlY_d+hlY_n(end) fliplr(lhY_d+lhY_n(end))], [156.5  206.5  194]./255,'EdgeColor','none','FaceAlpha',.3)
    
    patch([T_n fliplr(T_n)], [lhY_n fliplr(hhY_n)], [166 216 181]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([con_x fliplr(con_x)], [con_lhy fliplr(con_hhy)], [166 216 181]./255,'EdgeColor','none','FaceAlpha',.3)
    hold on
    patch([T_d+tnstop+twash fliplr(T_d+tnstop+twash)], [lhY_d+lhY_n(end) fliplr(hhY_d+hhY_n(end))], [166 216 181]./255,'EdgeColor','none','FaceAlpha',.3)
    
    
   %Error Bars
   hold on
   md= errorbar(x,y,yneg,ypos,'sb','MarkerSize',10,...
   'MarkerEdgeColor','black','MarkerFaceColor','black');
s
%Define Legend
lgd = legend([cphh cplh cphl cpll cp md],'Model 14.5 ug/mg, 30% D:W','Model 9 ug/mg, 30% D:W','Model 14.5 ug/mg, 20% D:W','Model 9 ug/mg, 20% D:W','Model Fitted Ym','Experimental Data','FontSize',16,'Location','Northeast','NumColumns',1);
hold off

%Save figures
saveas(pobj_ym(2,col-1),char(strcat('Figures/Ym/',c,'_2.pdf')));
saveas(pobj_ym(2,col-1),char(strcat('Figures/Ym/',c,'_2.fig')));
end