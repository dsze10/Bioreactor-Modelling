clc
clearvars
close all

%% PlateauSim with Test Parameters

%Plateau Sim code

AbsMax = 100;

TotalDNA = 100;
ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
Reactor_Volume = 25;

Abs2Dna = ConversionFactor*Reactor_Volume;
Dna2Abs = 1/Abs2Dna;

tn0 = 0;
tndur = 60;
tnend= tn0+tndur;
MaxN = 0.8;
kn= 0.01;
y0n= 0;
%ymn = TotalDNA*MaxN

twash = 20;


td0 = tnend+twash;
tddur = 60;
tdend = td0+tddur;
MaxD = 0.8;
kd= 0.05;

y0d = 0;
%ymd = (TotalDNA - yn(end))*MaxD

open('PlateauSim')
sim('PlateauSim')

%% Estimation of PlateauSim Parameters with Experimental Data

clearvars -except sheets Normdat Fitdat Mass NatDNA FinDNA  
Figures = {'Fig3.auVF1','Fig3.auVF2','Fig7.shVF.DS'};
Params = {'MaxN', 'kn', 'MaxD','kd'};
EstParam = array2table(zeros(length(Figures),length(Params)), 'VariableNames',Params,'RowNames',Figures);
%Parameters to estimate
        MaxN = 0.5;
        kn= 0.05;

        %ymn = TotalDNA*MaxN
        
        MaxD = 0.5;
        kd= 0.05;

        %ymd = (TotalDNA - yn(end))*MaxD
%         hws = get_param('PlateauSim_Opti', 'modelworkspace');
%         hws.assignin('MaxN', MaxN);
%         hws.assignin('MaxD', MaxD);
%         hws.assignin('kn', kn);
%         hws.assignin('kd', kd);
%Fixed Parameters
        y0n= 0;
        y0d = 0;
        D2W = 0.3; 
        DNA = 14.5;
        ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
        Reactor_Volume = 25;
        Abs2Dna = ConversionFactor*Reactor_Volume;
        Dna2Abs = 1/Abs2Dna;
        
%Get Data *need to run Excel_Import.m first
for c= Figures
        c= Figures{1};
        
        col = find(strcmp(Fitdat(1,:),c));
        NaDeoDat = Fitdat{5,col};
        NaDeoParam = Fitdat{6,col};
        DNaseDat = Fitdat{7,col};
        DNaseParam = Fitdat{8,col};

    %Constants
        Mass_ECM_wet = Mass{strcmp(Mass.Fig, c), 3}; %mg 
        Mass_ECM_dry = Mass_ECM_wet*D2W;
        TotalDNA =  DNA*Mass_ECM_dry; %#ok<NASGU>

    tn0 = 0;
        tndur= NaDeoParam{strcmp(NaDeoParam.Param, 'Dur'), 2}; %Duration of the Sodium Doexycholate Step in min 
        tnend= tn0+tndur; 
    twash = DNaseParam{strcmp(DNaseParam.Param, 'twash'), 2};
    td0 = tnend+twash;
        tddur=DNaseParam{strcmp(DNaseParam.Param, 'Dur'), 2}; %%Duration of the DNAse Step in min
        tdend = td0+tddur;

    %Format Raw Data
            ends = size(NaDeoDat);
            DNase = DNaseDat;
            DNase.('Time (min)')= DNase.('Time (min)')+twash+tnend;
            ScatDat = [NaDeoDat;DNase];
            ScatDat.Mean =  mean(ScatDat{:,3:ends(2)},2,'omitnan','double');
            ScatDat.Mean(ends(1)+1:end) = ScatDat.Mean(ends(1)+1:end)+ScatDat.Mean(ends(1));
          
AbsMax = (ceil(ScatDat.Mean(end))+1);
Ttotal = (ceil(tdend/10)+1)*10;

EstData = timeseries(ScatDat.Mean,ScatDat.('Time (min)'));
PlateauSim_Opti
set_param('PlateauSim_Opti', 'SimulationCommand', 'Update')
EstParam{c,:} = [MaxN, kn, MaxD,kd];

%Generated Optimization function
% [pOpt, Info] = parameterEstimation_PlateauSim_Opti(, EstData);

%After sim runs save Parameter data into table
EstParam{c,:} = [pOpt(2,1).Value, pOpt(4,1).Value, pOpt(1,1).Value, pOpt(3,1).Value];

%EstParam{:,strcat(c,'_Prism')} = {c NaDeoParam{strcmp(NaDeoParam.Param, 'YM'), 2} NaDeoParam{strcmp(NaDeoParam.Param, 'k'), 2} DNaseParam{strcmp(DNaseParam.Param, 'YM'), 2} DNaseParam{strcmp(DNaseParam.Param, 'k'), 2}}'
end

 writetable(EstParam,strcat('SimulinkEstParam.csv'),'WriteRowNames',true);
 %% Run with estimated Params
clearvars -except sheets Normdat Fitdat Mass NatDNA FinDNA EstParam
Figures = {'Fig3.auVF1','Fig3.auVF2','Fig7.shVF.DS'};
DNAvars = {'TotalDNA','DNAleft'};
DNAvals = array2table(zeros(length(Figures),length(DNAvars)), 'VariableNames',DNAvars,'RowNames',Figures);
 %Fixed Parameters
        y0n= 0;
        y0d = 0;
        D2W = 0.3; 
        DNA = 14.5;
        ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
        Reactor_Volume = 25;
        Abs2Dna = ConversionFactor*Reactor_Volume;
        Dna2Abs = 1/Abs2Dna;
 
  for c= Figures
        col = find(strcmp(Fitdat(1,:),c));
        NaDeoParam = Fitdat{6,col};
        DNaseParam = Fitdat{8,col};

    %Constants
        Mass_ECM_wet = Mass{strcmp(Mass.Fig, c), 3}; %mg 
        Mass_ECM_dry = Mass_ECM_wet*D2W;
        TotalDNA =  DNA*Mass_ECM_dry; %#ok<NASGU>

    tn0 = 0;
        tndur= NaDeoParam{strcmp(NaDeoParam.Param, 'Dur'), 2}; %Duration of the Sodium Doexycholate Step in min 
        tnend= tn0+tndur; 
    twash = DNaseParam{strcmp(DNaseParam.Param, 'twash'), 2};
    td0 = tnend+twash;
        tddur=DNaseParam{strcmp(DNaseParam.Param, 'Dur'), 2}; %%Duration of the DNAse Step in min
        tdend = td0+tddur;
        
    MaxN = EstParam{c,'MaxN'};
    kn= EstParam{c,'kn'};

    MaxD = EstParam{c,'MaxD'};
    kd= EstParam{c,'kd'};
    
   simOut = sim('PlateauSim_Opti');
        
  DNAvals{c,:} = [simOut.TotalDNA(end),simOut.DNAleft(end)];
  end
   ExpVals = [EstParam, DNAvals];
   writetable(ExpVals,strcat('SimulinkExpVals.csv'),'WriteRowNames',true);
        
 %% Optimization of Parameters
clearvars -except sheets Normdat Fitdat Mass NatDNA FinDNA EstParam
%Fixed Parameters
        y0n= 0;
        y0d = 0;
        D2W = 0.3; 
        DNA = 14.5;
        ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
        Reactor_Volume = 25;
        Abs2Dna = ConversionFactor*Reactor_Volume;
        Dna2Abs = 1/Abs2Dna;
  
EstParam = readtable(strcat('SimulinkEstParam.csv'));
%Set Extimated Parameters
cfit='Fig3.auVF2';
cdat = 'Fig3.auVF1';
    col= find(strcmp(Fitdat(1,:),cdat));
    NaDeoDat = Fitdat{5,col};
    NaDeoParam = Fitdat{6,col};
    DNaseDat = Fitdat{7,col};
    DNaseParam = Fitdat{8,col};

    MaxN = EstParam{cfit,'MaxN'};
    kn= EstParam{cfit,'kn'}{1,1};

    MaxD = EstParam{cfit,'MaxD'};
    kd= EstParam{cfit,'kd'};
    
%Set New Experimental Data
    Mass_ECM_wet = Mass{strcmp(Mass.Fig,cdat), 3}; %mg ;
    Mass_ECM_dry = Mass_ECM_wet*D2W;
    TotalDNA =  DNA*Mass_ECM_dry;
    
% Acceptable Min DNA
    MinDNA = 2*Mass_ECM_dry;

    tn0 = 0;
        tndur= NaDeoParam{strcmp(NaDeoParam.Param, 'Dur'), 2}; %Duration of the Sodium Doexycholate Step in min 
        tnend= tn0+tndur; 
    twash = DNaseParam{strcmp(DNaseParam.Param, 'twash'), 2};
    td0 = tnend+twash;
        tddur=DNaseParam{strcmp(DNaseParam.Param, 'Dur'), 2}; %%Duration of the DNAse Step in min
        tdend = td0+tddur;
        
sim('PlateauSim_Opti')
%% Extract Optimization Data

DVars = DesignVars.Data.Workspace.LocalWorkspace.DesignVars;
OptiDur = table({'Fig' 'tn0' 'tnend' 'td0' 'tdend'}');
OptiDur = table({'Fig' DVars(3,1).Name DVars(4,1).Name DVars(1,1).Name DVars(2,1).Name}');
OptiDur{:,strcat('Init_',cdat)} = {cdat tn0 tnend td0 tdend}';
OptiDur{:,strcat('Opti_',cdat)} = {cdat DVars(3,1).Value DVars(4,1).Value DVars(1,1).Value DVars(2,1).Value}';
writetable(OptiDur,strcat('SimulinkOptiDur','.csv'));

%% Response optimization 092522
clearvars -except sheets Normdat Fitdat Mass NatDNA FinDNA  
Figures = {'Fig3.auVF1','Fig3.auVF2','Fig7.shVF.DS'};
Params = {'MaxN', 'kn', 'MaxD','kd'};
Paramvals = {0.0921,0.0539,0.1392,0.0793;
             0.5888,0.0147,0.8037,0.0216;
             0.9966,0.0103,0.6375,0.0549};
EstParam = array2table(Paramvals, 'VariableNames',Params,'RowNames',Figures);

%Fixed Parameters
        y0n= 0;
        y0d = 0;
        D2W = 0.3; 
        DNA = 14.5;
        ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
        Reactor_Volume = 25;
        Abs2Dna = ConversionFactor*Reactor_Volume;
        Dna2Abs = 1/Abs2Dna;
        
cfit='Fig3.auVF2';
cfit='Fig7.shVF.DS';
cdat = 'Fig3.auVF2';
cdat ='Fig7.shVF.DS';
    col= find(strcmp(Fitdat(1,:),cdat));
    NaDeoDat = Fitdat{5,col};
    NaDeoParam = Fitdat{6,col};
    DNaseDat = Fitdat{7,col};
    DNaseParam = Fitdat{8,col};

    MaxN = EstParam{cfit,'MaxN'}{1,1};
    kn= EstParam{cfit,'kn'}{1,1};
    MaxD = EstParam{cfit,'MaxD'}{1,1};
    kd= EstParam{cfit,'kd'}{1,1};
    
%Set New Experimental Data
    Mass_ECM_wet = Mass{strcmp(Mass.Fig,cdat), 3}; %mg ;
    Mass_ECM_dry = Mass_ECM_wet*D2W;
    TotalDNA =  DNA*Mass_ECM_dry;

    tn0 = 0;
        tndur= NaDeoParam{strcmp(NaDeoParam.Param, 'Dur'), 2}; %Duration of the Sodium Doexycholate Step in min 
        tnend= tn0+tndur; 
    twash = DNaseParam{strcmp(DNaseParam.Param, 'twash'), 2};
    td0 = tnend+twash;
        tddur=DNaseParam{strcmp(DNaseParam.Param, 'Dur'), 2}; %%Duration of the DNAse Step in min
        tdend = td0+tddur;
        
% Acceptable Min DNA
    MinDNA = 2*Mass_ECM_dry;
    
sim('PlateauSim_Opti')

DVars = PlateauSim_Opti_sdosession.Data.Workspace.LocalWorkspace.Duration;
OptiDur = table({'Fig' 'tn0' 'tnend' 'td0' 'tdend'}');
OptiDur = table({'Fig' DVars(3,1).Name DVars(4,1).Name DVars(1,1).Name DVars(2,1).Name}');
OptiDur{:,strcat('Init_',cdat)} = {cdat tn0 tnend td0 tdend}';
OptiDur{:,strcat('Opti_',cdat)} = {cdat DVars(3,1).Value DVars(4,1).Value DVars(1,1).Value DVars(2,1).Value}';
writetable(OptiDur,strcat('SimulinkOptiDur','.csv'));
