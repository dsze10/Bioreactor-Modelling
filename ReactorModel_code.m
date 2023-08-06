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
            %writetable(EmpDat,filename,'Sheet',strcat(Normdat{1,c},'_',Normdat{r,1}));
            %writetable(EmpDat(:,2:end),strcat('NormCSV/',Normdat{1,c},'_',Normdat{r,1},'.csv'));
        end
    end
end

clearvars -except sheets Normdat


%% Import Fit from Prism 
%USING NORM DATA
%%%%
FitVal = readcell('Results042322.txt');

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
         %%%
        sheetname = find(strcmp(FitVal(:,1),strcat({'Nonlin fit of '},Fitdat{1,c},'_',Fitdat{n,1})));
         if ~isempty(sheetname)
             %%%
             Fitdat{n+1,c} = cell2table(FitVal(sheetname+2:sheetname+5,1:2),'VariableNames',{'Param' 'Val'});
             %%%
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

clearvars -except sheets Normdat Fitdat Mass NatDNA FinDNA  
%% Plateau Sim Code
%Fixed Parameters
        y0n= 0;
        y0d = 0;
        D2W = 0.3; 
        DNA = 14.5;
        ConversionFactor= 50; %Conversion factor from ABS 260 to [ug/ml]
        Reactor_Volume = 25;
        Abs2Dna = ConversionFactor*Reactor_Volume;
        Dna2Abs = 1/Abs2Dna;
        
%Estimate Parameters
Figures = {'Fig3.auVF1','Fig3.auVF2','Fig7.shVF.DS'};
Params = {'MaxN', 'kn', 'MaxD','kd'};


%Get Data for each Figure RUN ONCE PER FIGURE
c= 'Fig7.shVF.DS';
        
        col = find(strcmp(Fitdat(1,:),c));
        NaDeoDat = Fitdat{5,col};
        NaDeoParam = Fitdat{6,col};
        DNaseDat = Fitdat{7,col};
        DNaseParam = Fitdat{8,col};

    %Constants
        Mass_ECM_wet = Mass{strcmp(Mass.Fig, c), 3}; %mg 
        Mass_ECM_dry = Mass_ECM_wet*D2W;
        TotalDNA =  DNA*Mass_ECM_dry; 

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
%%Copy data into simulink for parameter estimation experiment
EstData = timeseries(ScatDat.Mean,ScatDat.('Time (min)'));

PlateauSim_Opti
set_param('PlateauSim_Opti', 'SimulationCommand', 'Update')

%Run simulation to plot data through model
sim('PlateauSim_Opti')

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

cdat = 'Fig3.auVF1';
    col= find(strcmp(Fitdat(1,:),cdat));
    NaDeoDat = Fitdat{5,col};
    NaDeoParam = Fitdat{6,col};
    DNaseDat = Fitdat{7,col};
    DNaseParam = Fitdat{8,col};

    
cfit='Fig3.auVF2';
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

