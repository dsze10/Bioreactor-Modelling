close all
clc
clear variables;
%% Import Excel file of Raw abs260 data
% Excel file format :Each sheet should contain one set of data at A1

% Stores each data from each sheet in column of cell array
%   - Column 1: ID (defines what data is extracted in the row)
% row 1: sheet names
% row 2: Full data
% row 3: searches for 'pbs' case sensitive (for fig 2 data)
% row 4: searches for 'NaDeo' case INsensitive
% row 5: searches for 'DNase' case INsensitive
% row 6: searches for 'DNase2' case INsensitive (for fig 7)

% gets sheet names from file 
                                  % change file name here %
sheets = ['ID',cellstr(sheetnames('102921_RawData.xlsx'))'];
 
for k=2:(numel(sheets))
  sheets{2,1} = 'RawTable';
  sheets{2,k} = readtable('102921_RawData.xlsx','sheet',k-1);
  Ind = find(contains(table2cell((sheets{2,k}(:,1))),'pbs'));
  sheets{3,1} = 'pbs';
  if ~isempty(Ind) && Ind(1) > 1
    sheets{3,k} = sheets{2,k}([Ind(1)-1;Ind],:);
  elseif ~isempty(Ind)
    sheets{3,k} = sheets{2,k}(Ind,:);
  else
    sheets{4,k} = [];
  end
  Ind = find(contains(table2cell((sheets{2,k}(:,1))),'NaDeo_','IgnoreCase',true));
  sheets{4,1} = 'NaDeo';
  if ~isempty(Ind) && Ind(1) > 1
      sheets{4,k} = sheets{2,k}([Ind(1)-1;Ind],:);
  elseif ~isempty(Ind) 
     sheets{4,k} = sheets{2,k}(Ind,:);
  else
     sheets{4,k} = [];
  end
  Ind = find(contains(table2cell((sheets{2,k}(:,1))),'DNase_','IgnoreCase',true));
  sheets{5,1} = 'DNase';
   if ~isempty(Ind) && Ind(1) > 1
       sheets{5,k} = sheets{2,k}([Ind(1)-1;Ind],:);
   elseif ~isempty(Ind) 
       sheets{5,k} = sheets{2,k}(Ind,:);
   else
       sheets{5,k} = [];
   end
  Ind = find(contains(table2cell((sheets{2,k}(:,1))),'DNase2_','IgnoreCase',true));
  sheets{6,1} = 'DNase2';
  if ~isempty(Ind) && Ind(1) > 1
      sheets{6,k} = sheets{2,k}([Ind(1)-1;Ind],:);
  elseif ~isempty(Ind) 
      sheets{6,k} = sheets{2,k}(Ind,:);
  else
      sheets{6,k} = [];
  end
end

clearvars -except sheets 
%% Normalization 
%in datatable reagent:time:Abs...
Normdat = sheets;
filename = 'NormData.xlsx';
for c=2:(size(Normdat,2))
    for r =3:(size(Normdat,1))
        if ~isempty(Normdat{r,c})
            Normdat{r,c}.Properties.VariableNames(2) = {'Time (min)'};
            %Find value of start time
            time = Normdat{r,c}{1,2};
            %find mean abs260 value at time 0
            norm = nanmean(Normdat{r,c}{1,3:end});
            %Remove norm value from all Reps  
            Normdat{r,c}{:,3:end} = (Normdat{r,c}{:,3:end}) - norm;
            %Set time to start at 0
            Normdat{r,c}{:,2} = (Normdat{r,c}{:,2}) - time;
            EmpDat = convertvars(Normdat{r,c},@isnumeric,@NanEmp);
            %Save as excel and csv file
            writetable(EmpDat,filename,'Sheet',strcat(Normdat{1,c},'_',Normdat{r,1}));
            writetable(EmpDat(:,2:end),strcat('NormCSV/',Normdat{1,c},'_',Normdat{r,1},'.csv'));
        end
    end
end

clearvars -except sheets Normdat

%% Import Fit from Prism 
%USING NORM DATA
FitVal = readcell('NormCSV/Results042322.txt');
Fitdat = Normdat;
m= [1000;1000;500;1000;500;1000;1000;1000];
nDNA= [11.1;7.1;11.1;11.1;7.1;7.1;11.1;11.1];
fDNA= [2.3;0.8;3.2;1.2;0.1;0.1;2.5;2.5];
Mass = table(Fitdat(1,2:end)',repmat({'Mass'},length(m),1),m,'VariableNames',{'Fig' 'Param' 'Val'});
NatDNA = table(Fitdat(1,2:end)',repmat({'NatDNA'},length(nDNA),1),nDNA,'VariableNames',{'Fig' 'Param' 'Val'});
FinDNA = table(Fitdat(1,2:end)',repmat({'FinDNA'},length(fDNA),1),fDNA,'VariableNames',{'Fig' 'Param' 'Val'});


FitTab = [{'Fig';'Reagent';'YM';'Y0';'k';'R squared'}];


%Insert empty rows
row = 3;
for r = 0:(size(Fitdat,1)-row) %not including PBS
    n=row+(r*2);
    Fitdat(n+1:end+1,:) = Fitdat(n:end,:);
    Fitdat(n+1,:) = {[]};
    Fitdat{n+1,1} = strcat(Fitdat{n,1},'_Param');
    for c=2:(size(Normdat,2))
        %find sheetname from Results.txt
         sheetname = find(strcmp(FitVal(:,1),strcat({'Nonlin fit of '},Fitdat{1,c},'_',Fitdat{n,1})));
         if ~isempty(sheetname)
             Fitdat{n+1,c} = cell2table(FitVal(sheetname+2:sheetname+5,1:2),'VariableNames',{'Param' 'Val'});
             
             FitTab = [FitTab,[Fitdat{1,c};Fitdat{n,1};FitVal(sheetname+2:sheetname+5,2)]];
             
             Dur = cell2table({'Dur', Fitdat{n,c}{end,2}},'VariableNames',{'Param' 'Val'});
             T0 = cell2table({'T0', sheets{row+r,c}{1,2}},'VariableNames',{'Param' 'Val'});
             Tend = cell2table({'Tend', sheets{row+r,c}{end,2}},'VariableNames',{'Param' 'Val'});
             Fitdat{n+1,c} = [Fitdat{n+1,c};Mass(c-1,2:3);NatDNA(c-1,2:3);FinDNA(c-1,2:3);T0;Tend;Dur]; %Mass unused
             
             if row+r >3 && ~isempty(sheets{row+r-1,c})&& ~isempty(sheets{row+r,c})
              Twash = cell2table({'twash', sheets{row+r,c}{1,2} - sheets{row+r-1,c}{end,2}},'VariableNames',{'Param' 'Val'});
              Fitdat{n+1,c} = [Fitdat{n+1,c};Twash];
             end
             if iscell(Fitdat{n+1,c}.Val)
                 Fitdat{n+1,c}.Val(cellfun(@ischar,Fitdat{n+1,c}.Val)) = {nan};
                 Fitdat{n+1,c}.Val = cell2mat(Fitdat{n+1,c}.Val); 
             end
         end  
    end
end

%writetable(FitTab,strcat('NormCSV/','FitVars.csv'));
%clearvars -except sheets Normdat Scaledat Fitdat Mass NatDNA FinDNA TMdat pobj

%% Run with Model with Fitted Ym
clc
close all

% Fixed Constants
t0=0;
ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
Reactor_Volume = 25; %ml
hD2W =0.3; 
lD2W =0.2; 
highDNA = 14.5;
lowDNA = 9;
nMax = 1;
dMax = 1;

pobj_ym = gobjects(2,size(Fitdat,2)-1);

%For all figures
for c= Fitdat(1,2:end)
    %{'Fig2.shVF','Fig3.auVF2','Fig7.shVF.DS'}
    %{'Fig2.shVF','Fig2.shSG','Fig3.auVF1','Fig3.auVF2','Fig3.auSG1','Fig3.auSG2','Fig7.shVF.DS'} 
    col = find(strcmp(Fitdat(1,:),c));
    
    NaDeoDat = Fitdat{5,col};
    NaDeoParam = Fitdat{6,col};
    DNaseDat = Fitdat{7,col};
    DNaseParam = Fitdat{8,col};
   
    %Get Mass
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
    T_n=t0:n_duration;
    
    llymn=(lltotal_DNA*nMax/(Reactor_Volume))/ConversionFactor; %Max absorbance
    llY_n=BioReactor_model(kn,y0n,llymn,T_n)*ConversionFactor* Reactor_Volume;
    
    hlymn=(hltotal_DNA*nMax/(Reactor_Volume))/ConversionFactor; %Max absorbance
    hlY_n=BioReactor_model(kn,y0n,hlymn,T_n)*ConversionFactor* Reactor_Volume;
    
    lhymn=(lhtotal_DNA*nMax/(Reactor_Volume))/ConversionFactor; %Max absorbance
    lhY_n=BioReactor_model(kn,y0n,lhymn,T_n)*ConversionFactor* Reactor_Volume;
    
    hhymn=(hhtotal_DNA*nMax/(Reactor_Volume))/ConversionFactor; %Max absorbance
    hhY_n=BioReactor_model(kn,y0n,hhymn,T_n)*ConversionFactor* Reactor_Volume;
    
    ymn=NaDeoParam{strcmp(NaDeoParam.Param, 'YM'), 2};
    Y_n=BioReactor_model(kn,y0n,ymn,T_n)*ConversionFactor* Reactor_Volume;
    
%Get Twash
    twash = DNaseParam{strcmp(DNaseParam.Param, 'twash'), 2};

%For DNAse
    tdstop=tnstop+twash+d_duration; %min
    kd= DNaseParam{strcmp(DNaseParam.Param, 'k'), 2}; %Overall rate constant min^-1
    y0d=0; % Initial DNA concentration [ug/ml] after PBS wash
    T_d=t0:d_duration;
        
    llDNA_left=lltotal_DNA*dMax-llY_n(end);
    llymd=(llDNA_left/Reactor_Volume)/ConversionFactor;
    llY_d=BioReactor_model(kd,y0d,llymd,T_d)*ConversionFactor*Reactor_Volume;
    
    hlDNA_left=hltotal_DNA*dMax-hlY_n(end);
    hlymd=(hlDNA_left/Reactor_Volume)/ConversionFactor;
    hlY_d=BioReactor_model(kd,y0d,hlymd,T_d)*ConversionFactor*Reactor_Volume;
    
    lhDNA_left=lhtotal_DNA*dMax-lhY_n(end);
    lhymd=(lhDNA_left/Reactor_Volume)/ConversionFactor;
    lhY_d=BioReactor_model(kd,y0d,lhymd,T_d)*ConversionFactor*Reactor_Volume;
    
    hhDNA_left=hhtotal_DNA*dMax-hhY_n(end);
    hhymd=(hhDNA_left/Reactor_Volume)/ConversionFactor;
    hhY_d=BioReactor_model(kd,y0d,hhymd,T_d)*ConversionFactor*Reactor_Volume;
    
    ymd=DNaseParam{strcmp(DNaseParam.Param, 'YM'), 2};
    Y_d=BioReactor_model(kd,y0d,ymd,T_d)*ConversionFactor*Reactor_Volume;
    
 % Get Scatterplot Data for Plot 1

    cend = size(NaDeoDat,2);
    DNase = DNaseDat;
    DNase.('Time (min)')= DNase.('Time (min)')+twash+tnstop;
    ScatDat = [NaDeoDat;DNase];
    ScatDat.Mean =  mean(ScatDat{:,3:cend},2,'omitnan','double');
    ScatDat.STD = std(ScatDat{:,3:cend},0,2,'omitnan');
    ScatDat.SEM = ScatDat.STD./sqrt(sum(~isnan(ScatDat{:,3:cend}),2));
   
    
    
%Plot Figure 1
    x = ScatDat.('Time (min)')';
    y = ScatDat.Mean'*ConversionFactor* Reactor_Volume;
    yneg = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
    ypos = ScatDat.SEM'*ConversionFactor* Reactor_Volume;

    
    pobj_ym(1,col-1) = figure('visible','on');
    set(gca,'FontSize',16);
    ylabel('DNA [ug]','FontSize',20);
    xlabel('Time (min)','FontSize',20);
    title(strcat(c,' (',num2str(Mass_ECM_wet),' mg)'),'FontSize',24);
    ylim([0 inf]);
    xlim([0 inf]);
    
    
    con_x = [T_n(end) T_n(end)+twash];
    

    hold on
    plot(T_n,llY_n,...
        T_d+tnstop+twash,llY_d,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);
    hold on
    con_lly = [llY_n(end) llY_d(1)];
    cpll= plot(con_x,con_lly,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);

    
    hold on
    plot(T_n,hlY_n,...
        T_d+tnstop+twash,hlY_d,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    hold on
    con_hly = [hlY_n(end) hlY_d(1)];
    cphl= plot(con_x,con_hly,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    
    hold on
    plot(T_n,lhY_n,...
        T_d+tnstop+twash,lhY_d,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    hold on
    con_lhy = [lhY_n(end) lhY_d(1)];
    cplh= plot(con_x,con_lhy,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);

    
    hold on
    plot(T_n,hhY_n,...
        T_d+tnstop+twash,hhY_d,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    hold on
    con_hhy = [hhY_n(end) hhY_d(1)];
    cphh= plot(con_x,con_hhy,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    
    hold on
    plot(T_n,Y_n,...
        T_d+tnstop+twash,Y_d,...
        'LineWidth',3,'LineStyle','--','Color','b');
    hold on
    con_y = [Y_n(end) Y_d(1)];
    cp= plot(con_x,con_y,...
        'LineWidth',3,'LineStyle','--','Color','b');

    hold on

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
    
lgd = legend([cphh cplh cphl cpll cp md],'Model 14.5 ug/mg, 30% D:W','Model 9 ug/mg, 30% D:W','Model 14.5 ug/mg, 20% D:W','Model 9 ug/mg, 20% D:W','Model Fitted Ym','Measured Data','FontSize',16,'Location','Northeast','NumColumns',1);
hold off

saveas(pobj_ym(1,col-1),char(strcat('Figures/Ym/',c,'_1.pdf')));
%saveas(pobj_ym(1,col-1),char(strcat('Figures/Ym/',c,'_1.fig')));


% Get Scatterplot Data for Plot 2

    cend = size(NaDeoDat,2);
    NaDeoend = mean(NaDeoDat{end,3:cend},'omitnan','double');
    DNase = DNaseDat;
    DNase.('Time (min)')= DNase.('Time (min)')+twash+tnstop;
    DNase{:,3:cend} = DNase{:,3:cend}+NaDeoend;
    ScatDat = [NaDeoDat;DNase];
    ScatDat.Mean =  mean(ScatDat{:,3:cend},2,'omitnan','double');
    ScatDat.STD = std(ScatDat{:,3:cend},0,2,'omitnan');
    ScatDat.SEM = ScatDat.STD./sqrt(sum(~isnan(ScatDat{:,3:cend}),2));

%Plot Figure 2
    x = ScatDat.('Time (min)')';
    y = ScatDat.Mean'*ConversionFactor* Reactor_Volume;
    yneg = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
    ypos = ScatDat.SEM'*ConversionFactor* Reactor_Volume;

    con_x = [T_n(end) T_n(end)+twash];
    
    pobj_ym(2,col-1) = figure('visible','on');
    %Set Graphing param
    set(gca,'FontSize',16);
    ylabel('DNA [ug]','FontSize',20);
    xlabel('Time (min)','FontSize',20);
    title(strcat(c,' (',num2str(Mass_ECM_wet),' mg)'),'FontSize',24);
    ylim([0 inf]);
    xlim([0 inf]);
    
    hold on
    plot(T_n,llY_n,...
        T_d+tnstop+twash,llY_d+llY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);
    hold on
    con_lly = [llY_n(end) llY_n(end)];
    cpll= plot(con_x,con_lly,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);

    
    hold on
    plot(T_n,hlY_n,...
        T_d+tnstop+twash,hlY_d+hlY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    hold on
    con_hly = [hlY_n(end) hlY_n(end)];
    cphl= plot(con_x,con_hly,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    
    hold on
    plot(T_n,lhY_n,...
        T_d+tnstop+twash,lhY_d+lhY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    hold on
    con_lhy = [lhY_n(end) lhY_n(end)];
    cplh= plot(con_x,con_lhy,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    
    hold on
    plot(T_n,hhY_n,...
        T_d+tnstop+twash,hhY_d+hhY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    hold on
    con_hhy = [hhY_n(end) hhY_n(end)];
    cphh= plot(con_x,con_hhy,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    
    hold on
    plot(T_n,Y_n,...
        T_d+tnstop+twash,Y_d+Y_n(end),...
        'LineWidth',3,'LineStyle','--','Color','b');
    hold on
    con_y = [Y_n(end) Y_n(end)];
    cp= plot(con_x,con_y,...
        'LineWidth',3,'LineStyle','--','Color','b');

    hold on

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
    
lgd = legend([cphh cplh cphl cpll cp md],'Model 14.5 ug/mg, 30% D:W','Model 9 ug/mg, 30% D:W','Model 14.5 ug/mg, 20% D:W','Model 9 ug/mg, 20% D:W','Model Fitted Ym','Measured Data','FontSize',16,'Location','Northeast','NumColumns',1);
hold off

saveas(pobj_ym(2,col-1),char(strcat('Figures/Ym/',c,'_2.pdf')));
%saveas(pobj_ym(2,col-1),char(strcat('Figures/Ym/',c,'_2.fig')));
end
%set(pobj_ym(:,:),'Visible','off')
clearvars -except sheets Normdat Scaledat Fitdat Mass NatDNA FinDNA TMdat pobj pobj_ym
%% Testing Fig3 1000mg fit with 500mg
%plot
%For all figur
    
clc
close all
% Fixed Constants
t0=0;
ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
Reactor_Volume = 25; %ml
hD2W =0.3; 
lD2W =0.2; 
highDNA = 14.5;
lowDNA = 9;

Figures = {'Fig3.auVF2','Fig3.auSG2','Fig7.shVF.DS','Fig3.auVF2'};
datFig = {'Fig3.auVF1','Fig3.auSG1','Fig3.auVF2','Fig7.shVF.DS'};
    %{'Fig2.shVF','Fig2.shSG','Fig3.auVF1','Fig3.auVF2','Fig3.auSG1','Fig3.auSG2','Fig7.shVF','Fig7.shVF.DS'}
Names = erase(append(Figures,'_',datFig),"Fig");
rmDNArows = ["kn", "NADeo YM", "kd", "DNase YM","rmDNA Pico Per mg","rmDNA Pico, 30% D:W", "rmDNA Pico, 20% D:W","rmDNA rawData", "rmDNA Ym", "rmDNA 14.5 ug/mg, 30% D:W", "rmDNA 9 ug/mg, 30% D:W", "rmDNA 14.5 ug/mg, 20% D:W", "rmDNA 9 ug/mg, 20% D:W"];
rmDNAdat = array2table(zeros(length(rmDNArows),length(Figures)), 'VariableNames',Names,'RowNames',rmDNArows);

pobj_rm = gobjects(2,size(Names,2));
%plot
for cn= 1:length(Figures)
    c=Figures(cn);
    col = find(strcmp(Fitdat(1,:),c));
    
    cdat = datFig(cn);
    coldat = find(strcmp(Fitdat(1,:),cdat));
    
    NaDeoDat = Fitdat{5,coldat};
    NaDeoParam = Fitdat{6,col};
    DNaseDat = Fitdat{7,coldat};
    DNaseParam = Fitdat{8,col};
   
    %Get Mass
    Mass_ECM_wet =  Mass{strcmp(Mass.Fig, cdat), 3}; %mg
    lMass_ECM_dry = Mass_ECM_wet*lD2W; 
    hMass_ECM_dry = Mass_ECM_wet*hD2W; 
    lltotal_DNA =  lowDNA*lMass_ECM_dry;
    hltotal_DNA =  highDNA*lMass_ECM_dry;
    lhtotal_DNA =  lowDNA*hMass_ECM_dry;
    hhtotal_DNA =  highDNA*hMass_ECM_dry; 

 %Get plotting Data
    %FitDat row 6:NaDeo ; row 8: Dnase
    n_duration= Fitdat{6,coldat}{strcmp(Fitdat{6,coldat}.Param, 'Dur'), 2}; %Duration of the Sodium Doexycholate Step in min 
    d_duration=Fitdat{8,coldat}{strcmp(Fitdat{8,coldat}.Param, 'Dur'), 2}; %%Duration of the DNAse Step in min
       
%For NaDeoxy
    tnstop=t0+n_duration; %min
    kn= Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'k'), 2}; %Overall rate constant min^-1
    y0n=0; % Initial DNA concentration [ug/ml] after PBS wash
    T_n=t0:n_duration;
    
    llymn=(lltotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    llY_n=BioReactor_model(kn,y0n,llymn,T_n)*ConversionFactor* Reactor_Volume;
    
    hlymn=(hltotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    hlY_n=BioReactor_model(kn,y0n,hlymn,T_n)*ConversionFactor* Reactor_Volume;
    
    lhymn=(lhtotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    lhY_n=BioReactor_model(kn,y0n,lhymn,T_n)*ConversionFactor* Reactor_Volume;
    
    hhymn=(hhtotal_DNA/(Reactor_Volume))/ConversionFactor; %Max absorbance
    hhY_n=BioReactor_model(kn,y0n,hhymn,T_n)*ConversionFactor* Reactor_Volume;
    
    ymn=Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'YM'), 2}/(Mass{strcmp(Mass.Fig, c), 3}/Mass{strcmp(Mass.Fig, cdat), 3});
    Y_n=BioReactor_model(kn,y0n,ymn,T_n)*ConversionFactor* Reactor_Volume;
    
%Get Twash
    twash = Fitdat{8,coldat}{strcmp(Fitdat{8,coldat}.Param, 'twash'), 2};

%For DNAse
    tdstop=tnstop+twash+d_duration; %min
    kd= Fitdat{8,col}{strcmp(Fitdat{8,col}.Param, 'k'), 2}; %Overall rate constant min^-1
    y0d=0; % Initial DNA concentration [ug/ml] after PBS wash
    T_d=t0:d_duration;
        
    llDNA_left=lltotal_DNA-llY_n(end);
    llymd=(llDNA_left/Reactor_Volume)/ConversionFactor;
    llY_d=BioReactor_model(kd,y0d,llymd,T_d)*ConversionFactor*Reactor_Volume;
    
    hlDNA_left=hltotal_DNA-hlY_n(end);
    hlymd=(hlDNA_left/Reactor_Volume)/ConversionFactor;
    hlY_d=BioReactor_model(kd,y0d,hlymd,T_d)*ConversionFactor*Reactor_Volume;
    
    lhDNA_left=lhtotal_DNA-lhY_n(end);
    lhymd=(lhDNA_left/Reactor_Volume)/ConversionFactor;
    lhY_d=BioReactor_model(kd,y0d,lhymd,T_d)*ConversionFactor*Reactor_Volume;
    
    hhDNA_left=hhtotal_DNA-hhY_n(end);
    hhymd=(hhDNA_left/Reactor_Volume)/ConversionFactor;
    hhY_d=BioReactor_model(kd,y0d,hhymd,T_d)*ConversionFactor*Reactor_Volume;
    
    ymd=Fitdat{8,col}{strcmp(Fitdat{8,col}.Param, 'YM'), 2}/(Mass{strcmp(Mass.Fig, c), 3}/Mass{strcmp(Mass.Fig, cdat), 3});
    Y_d=BioReactor_model(kd,y0d,ymd,T_d)*ConversionFactor*Reactor_Volume;

 % Get Scatterplot Data for Plot 1

    cend = size(NaDeoDat,2);
    DNase = DNaseDat;
    DNase.('Time (min)')= DNase.('Time (min)')+twash+tnstop;
    ScatDat = [NaDeoDat;DNase];
    ScatDat.Mean =  mean(ScatDat{:,3:cend},2,'omitnan','double');
    ScatDat.STD = std(ScatDat{:,3:cend},0,2,'omitnan');
    ScatDat.SEM = ScatDat.STD./sqrt(sum(~isnan(ScatDat{:,3:cend}),2));
   
    
    
%Plot Figure 1
    x = ScatDat.('Time (min)')';
    y = ScatDat.Mean'*ConversionFactor* Reactor_Volume;
    yneg = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
    ypos = ScatDat.SEM'*ConversionFactor* Reactor_Volume;

    
    pobj_rm(1,cn) = figure('visible','on');
    set(gca,'FontSize',16);
    ylabel('DNA [ug]','FontSize',20);
    xlabel('Time (min)','FontSize',20);
    title(strcat(c," Fit - ",cdat," data"),'FontSize',24);
    ylim([0 inf]);
    xlim([0 inf]);
    
    
    con_x = [T_n(end) T_n(end)+twash];
    

    hold on
    plot(T_n,llY_n,...
        T_d+tnstop+twash,llY_d,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);
    hold on
    con_lly = [llY_n(end) llY_d(1)];
    cpll= plot(con_x,con_lly,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);

    
    hold on
    plot(T_n,hlY_n,...
        T_d+tnstop+twash,hlY_d,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    hold on
    con_hly = [hlY_n(end) hlY_d(1)];
    cphl= plot(con_x,con_hly,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    
    hold on
    plot(T_n,lhY_n,...
        T_d+tnstop+twash,lhY_d,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    hold on
    con_lhy = [lhY_n(end) lhY_d(1)];
    cplh= plot(con_x,con_lhy,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);

    
    hold on
    plot(T_n,hhY_n,...
        T_d+tnstop+twash,hhY_d,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    hold on
    con_hhy = [hhY_n(end) hhY_d(1)];
    cphh= plot(con_x,con_hhy,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    
    hold on
    plot(T_n,Y_n,...
        T_d+tnstop+twash,Y_d,...
        'LineWidth',3,'LineStyle','--','Color','b');
    hold on
    con_y = [Y_n(end) Y_d(1)];
    cp= plot(con_x,con_y,...
        'LineWidth',3,'LineStyle','--','Color','b');

    hold on

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
   md= errorbar(x,y,yneg,ypos,'-sb','MarkerSize',10,...
   'MarkerEdgeColor','blue','MarkerFaceColor','blue');
    
lgd = legend([cphh cplh cphl cpll cp md],'Model 14.5 ug/mg, 30% D:W','Model 9 ug/mg, 30% D:W','Model 14.5 ug/mg, 20% D:W','Model 9 ug/mg, 20% D:W','Model Ym','Measured Data','FontSize',16,'Location','southoutside','NumColumns',1);
hold off

saveas(pobj_rm(1,cn),char(strcat('Figures/rmDNA/',Names(cn),'_1.pdf')));
%saveas(pobj_rm(1,cn),char(strcat('Figures/rmDNA/',Names(cn),'_1.fig')));

    
    
% Get Scatterplot Data for Plot 2

    cend = size(NaDeoDat,2);
    NaDeoend = mean(NaDeoDat{end,3:cend},'omitnan','double');
    DNase = DNaseDat;
    DNase.('Time (min)')= DNase.('Time (min)')+twash+tnstop;
    DNase{:,3:cend} = DNase{:,3:cend}+NaDeoend;
    ScatDat = [NaDeoDat;DNase];
    ScatDat.Mean =  mean(ScatDat{:,3:cend},2,'omitnan','double');
    ScatDat.STD = std(ScatDat{:,3:cend},0,2,'omitnan');
    ScatDat.SEM = ScatDat.STD./sqrt(sum(~isnan(ScatDat{:,3:cend}),2));

%Plot Figure 2
    x = ScatDat.('Time (min)')';
    y = ScatDat.Mean'*ConversionFactor* Reactor_Volume;
    yneg = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
    ypos = ScatDat.SEM'*ConversionFactor* Reactor_Volume;

    con_x = [T_n(end) T_n(end)+twash];
    
    pobj_rm(2,cn) = figure('visible','on');
    %Set Graphing param
    set(gca,'FontSize',16);
    ylabel('DNA [ug]','FontSize',20);
    xlabel('Time (min)','FontSize',20);
    title(strcat(c," fit - ",cdat," data"),'FontSize',24);
    ylim([0 inf]);
    xlim([0 inf]);
    
    hold on
    plot(T_n,llY_n,...
        T_d+tnstop+twash,llY_d+llY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);
    hold on
    con_lly = [llY_n(end) llY_n(end)];
    cpll= plot(con_x,con_lly,...
        'LineWidth',3,'LineStyle','-','Color',[140 190 215]./255);

    
    hold on
    plot(T_n,hlY_n,...
        T_d+tnstop+twash,hlY_d+hlY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    hold on
    con_hly = [hlY_n(end) hlY_n(end)];
    cphl= plot(con_x,con_hly,...
        'LineWidth',3,'LineStyle','-','Color',[150 200 202.5]./255);
    
    hold on
    plot(T_n,lhY_n,...
        T_d+tnstop+twash,lhY_d+lhY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    hold on
    con_lhy = [lhY_n(end) lhY_n(end)];
    cplh= plot(con_x,con_lhy,...
        'LineWidth',3,'LineStyle','-','Color',[163 213 185.5]./255);
    
    hold on
    plot(T_n,hhY_n,...
        T_d+tnstop+twash,hhY_d+hhY_n(end),...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    hold on
    con_hhy = [hhY_n(end) hhY_n(end)];
    cphh= plot(con_x,con_hhy,...
        'LineWidth',3,'LineStyle','-','Color',[173 223 173]./255);
    
        hold on
    plot(T_n,Y_n,...
        T_d+tnstop+twash,Y_d+Y_n(end),...
        'LineWidth',3,'LineStyle','--','Color','b');
    hold on
    con_y = [Y_n(end) Y_n(end)];
    cp= plot(con_x,con_y,...
        'LineWidth',3,'LineStyle','--','Color','b');

    hold on

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
   md= errorbar(x,y,yneg,ypos,'-sb','MarkerSize',10,...
   'MarkerEdgeColor','blue','MarkerFaceColor','blue');

   lgd = legend([cphh cplh cphl cpll cp md],'Model 14.5 ug/mg, 30% D:W','Model 9 ug/mg, 30% D:W','Model 14.5 ug/mg, 20% D:W','Model 9 ug/mg, 20% D:W','Model Ym','Measured Data','FontSize',16,'Location','southoutside','NumColumns',1);
 

hold off

saveas(pobj_rm(2,cn),char(strcat('Figures/rmDNA/',Names(cn),'_2.pdf')));
%saveas(pobj_rm(2,cn),char(strcat('Figures/rmDNA/',Names(cn),'_2.fig')));

    dpdiff = Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'NatDNA'), 2}-Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'FinDNA'), 2};
    ltpdiff = (Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'NatDNA'), 2}-Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'FinDNA'), 2})*lMass_ECM_dry;
    htpdiff = (Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'NatDNA'), 2}-Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'FinDNA'), 2})*hMass_ECM_dry;
   
    rdiff = y(end);
    mdiff = (Y_d(end)+Y_n(end));
    lldiff = (llY_d(end)+llY_n(end));
    lhdiff = (lhY_d(end)+lhY_n(end));
    hldiff = (hlY_d(end)+hlY_n(end));
    hhdiff = (hhY_d(end)+hhY_n(end));


    rmDNAcol = cell2mat({kn; Fitdat{6,col}{strcmp(Fitdat{6,col}.Param, 'YM'), 2}; ...
                kd; Fitdat{8,col}{strcmp(Fitdat{8,col}.Param, 'YM'), 2}; ...
                dpdiff; htpdiff; ltpdiff; rdiff; ...
                mdiff; hhdiff;  hldiff; lhdiff;  lldiff});
       
    rmDNAdat{:,Names(cn)} = rmDNAcol;
end
%set(pobj_rm(:,:),'Visible','off')
writetable(rmDNAdat,'rmDNA_Data.csv','WriteRowNames',true);

clearvars -except sheets Normdat Scaledat Fitdat Mass NatDNA FinDNA TMdat pobj pobj_ym pobj_rm rmDNAdat

%% lsqcurvefit model in Matlab
clc
close all
% Fixed Constants
t0=0;
ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
Reactor_Volume = 25; %ml

%Set Figure and Data
    col = 9;
    cend = size(Fitdat{5,col},2); %Last column of Abs data
    
    NaDeoDat = Fitdat{5,col};
    NaDeoParam = Fitdat{6,col};
    DNaseDat = Fitdat{7,col};
    DNaseParam = Fitdat{8,col};
    
%Scatterplot data for NaDeo
    ScatDat = NaDeoDat;
    ScatDat.Mean =  mean(ScatDat{:,3:cend},2,'omitnan','double');
    ScatDat.STD = std(ScatDat{:,3:cend},0,2,'omitnan');
    ScatDat.SEM = ScatDat.STD./sqrt(sum(~isnan(ScatDat{:,3:cend}),2));

    xscat = ScatDat.('Time (min)')';
    yscat = ScatDat.Mean'*ConversionFactor* Reactor_Volume;
    yneg = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
    ypos = ScatDat.SEM'*ConversionFactor* Reactor_Volume;
        
%Get Fitting Data
    Stackdata =  stack(NaDeoDat,3:cend,'NewDataVariableName','Abs',...
          'IndexVariableName','Rep');
    Stackdata = rmmissing(Stackdata);
    tdata = Stackdata{:,2};
    ydata = Stackdata{:,4};
    
   %Plateau
   fun = @(t,tdata)t(1)-((t(2)-t(1)).*exp((-t(3)).*tdata));
   
   x0 = [0,0,0];
   x= lsqcurvefit(fun,x0,tdata,ydata);

   Fit = figure();
   plot(tdata,fun(x,tdata)*ConversionFactor* Reactor_Volume,'LineWidth',2)
   hold on
   hold on
   errorbar(xscat,yscat,yneg,ypos,'-s','MarkerSize',5,...
   'MarkerEdgeColor',[ 0.5843 0.8157 0.9882],'MarkerFaceColor',[ 0.5843 0.8157 0.9882], 'LineWidth',1)
lgd = legend('Lsqfit','Measured Data','FontSize',16,'Location','southoutside');