function [Y] = DATA_COVID_UK_Prelims(plot_prelims)
 
% Data retrieval function for COVID modelling
% Using PHE data
% https://coronavirus.data.gov.uk/
%
% with census data from each region e.g.
% https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/censusoutputareaestimatesintheyorkshireandthehumberregionofengland
%
% and transform lower area code to local authority code
% https://data.gov.uk/dataset/ec39697d-e7f4-4419-a146-0b9c9c15ee06/output-area-to-lsoa-to-msoa-to-local-authority-district-december-2017-lookup-with-area-classifications-in-great-britain
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging
% FROM
% Karl Friston
% $Id: DATA_COVID_US.m 7838 2020-04-23 17:40:45Z karl $
 

% load population meta data
load 'MetadataUK.mat'
load 'ages_per_LA_table'  % source 2018 census
load 'dwellings_regional' %  Source: 2019 survery - Office for National Statistics - Labour Force Survey (LFS)
%--------------------------------------------------------------------------

% load data from https://coronavirus.data.gov.uk/
%--------------------------------------------------------------------------

try
    C  = importdata('coronavirus-cases_latest.csv',',',1);
    D  = importdata('coronavirus-deaths_latest.csv',',',1); % note England figures are called UK deaths from 15/03 to 27/03 - adjust below - England and 122 deaths on March 26th
  
catch
    
    clc, warning('Please load csv files (census and cases) into the current working directory/paths set')
    help DATA_COVID_UK 
end

 

nations   = find(strcmp(C.textdata(:,3),'Nation'));
uniqueN   = unique(C.textdata(nations,1));
 
regions   = find(strcmp(C.textdata(:,3),'Region'));
uniqueR   = unique(C.textdata(regions,1));
 
LAs       = find(strcmp(C.textdata(:,3),'Upper tier local authority'));
uniLAname = unique(C.textdata(LAs,1));


% Organise days - adjust for  May with months_days([2:4] etc. First_date : {'2020-01-30'};
%--------------------------------------------------------------------------

% cases
% months        =  {'01'  ,'02' ,  '03' , '04' ,'05' };                   % -- Jan -- Feb -- March -- April -- May
% month_days    =  [30:31  1:29   1:31   1:30  1:18];                     % -- Jan -- Feb -- March -- April -- May
% 
% day_beginning =  [30   1  1    1   1];
% day_ending    =  [31  29  31   30  18];

months        =  {'02' ,  '03' , '04' ,'05' };                   % -- Jan -- Feb -- March -- April -- May
month_days    =  [ 10:29   1:31   1:30  1:18];                     % -- Jan -- Feb -- March -- April -- May

day_beginning =  [ 10  1    1   1];
day_ending    =  [ 29  31   30  18];

i = 1;
for m = 1:length(months)
    for d = day_beginning(m):day_ending(m)
        if month_days(i) < 10
            all_dates{i} =['2020-',months{m},'-0',num2str(month_days(i))];
        else
            all_dates{i} =['2020-',months{m},'-',num2str(month_days(i))];
        end
        i= i+1;
    end
end
dates_recorded = unique(C.textdata(2:end,4) ) ;

% Deaths in England from {'2020-03-27'} day 
first_death_recordings = find(strcmp('2020-03-07',all_dates))

% Get National Time Series  ----  Deaths and Cases
%--------------------------------------------------------------------------
Y.Nation.name   = uniqueN;
Y.Nation.deaths = zeros(length(all_dates),1);
Y.Nation.cases  = zeros(length(all_dates),1);
Y.Nation.pop    = sum(RegionT.pop)/1e6;              % population in millions

for i = 1:length(all_dates)
    
    day_c_all           = find(strcmp(all_dates{i},C.textdata(:,4)));
    day_c_rel           = find(strcmp(Y.Nation.name,C.textdata(day_c_all,1)));
    if day_c_rel
        if isnan(C.data(day_c_all(day_c_rel)-1,1)) % some missing entries in csv file
              Y.Nation.cases(i)   = 0;
        else
             Y.Nation.cases(i)   = C.data(day_c_all(day_c_rel)-1,1);
        end
    end
    
    if i > first_death_recordings
        day_d_all           = find(strcmp(all_dates{i},D.textdata(:,4)));
        day_d_rel           = find(strcmp(Y.Nation.name,D.textdata(day_d_all,1)));
        if day_d_rel
              if isnan(D.data(day_d_all(day_d_rel)-1,1))
                  Y.Nation.deaths(i) = 0;
              else
                    Y.Nation.deaths(i)  = D.data(day_d_all(day_d_rel)-1,1);
              end
        end
    end
    
    
end

 code_ages = find(agesperLA.Code  == 92000001); % England
 for ages = 1:91
            Y.Nation.ages_0_90(ages) = eval(['agesperLA.VarName',num2str(ages+4),'(',num2str(code_ages),')']);
 end
        
 
% Get Regional and Subregional Time Series ---- Cases Only
%--------------------------------------------------------------------------
region_ordered_for_UTLA = {'London', 'Yorkshire_The_Humber','East_Midlands',  'West_Midlands' , 'East',...
    'South_West', 'South_East', 'North_West', 'North_East'}      ;
    
region_dwelling_table = {'London', 'YorkshireandTheHumber','EastMidlands',  'WestMidlands' , 'East',...
    'SouthWest', 'SouthEast', 'NorthWest', 'NorthEast'}      ;
 
  ages_Bournemouth_Christhchurh_Poole = [395784	3868	4087	4280	4267	4361	4438	4774	4545	4391	4396	4490	4185	4028	3920	3954	3769	3602	3945	4216	5770	6416	6092	5447	4571	4628	4501	4779	4765	4516	4570	4745	4939	4873	5179	4995	4987	5299	5644	5441	5152	4658	4599	4631	4820	4803	4964	5254	5368	5295	5236	5310	5569	5527	5508	5500	5327	5276	4951	4713	4641	4547	4550	4252	4114	4301	4157	4165	4325	4459	4442	4894	5470	4080	3971	3869	3583	3168	2680	2830	2853	2844	2646	2423	2350	2161	2034	1867	1785	1600	1315	5274];
  
for Reg = 1:9
    
    Y.Region(Reg).name    = RegionT.regions(Reg);
    Y.Region(Reg).cases   = zeros(length(all_dates),1);
    Y.Region(Reg).pop     = RegionT.pop(Reg)/1e6;         % population in millions
    
    
    for i = 1:length(all_dates)
        day_c_all                = find(strcmp(all_dates{i},C.textdata(:,4)));
        day_c_rel                = find(strcmp(Y.Region(Reg).name,C.textdata(day_c_all,1)));
        if day_c_rel
             if isnan(C.data(day_c_all(day_c_rel)-1,1)) % some missing entries in csv file
                  Y.Region(Reg).cases(i)   = 0;
             else
                  Y.Region(Reg).cases(i)   = C.data(day_c_all(day_c_rel)-1,1);
             end
        end
    end
    
    % --- subregion cases only
    subregi =  find(strcmp(region_ordered_for_UTLA{Reg} ,cellstr(RegionUTLAMetaData.Region)));
    for s = 1:length(subregi)
        
        Y.Region(Reg).subregion(s).cases = zeros(length(all_dates),1);
        Y.Region(Reg).subregion(s).pop   = RegionUTLAMetaData.Population(subregi(s));
        Y.Region(Reg).subregion(s).name  = cellstr(RegionUTLAMetaData.AreaUTLA(subregi(s)));
        Y.Region(Reg).subregion(s).code  = RegionUTLAMetaData.CodeUTLA(subregi(s));
       
        Table_code = Y.Region(Reg).subregion(s).code ;%10000000
    
        
        if  Table_code  <10000000 % add E0 to codes
            Adj_Table_code =['E0',num2str(Table_code)];
        else
            Adj_Table_code =['E',num2str(Table_code)];
        end
        
        for i = 1:length(all_dates)
            day_c_all                             = find(strcmp(all_dates{i},C.textdata(:,4)));
            day_c_rel                             = find(strcmp( Adj_Table_code,C.textdata(day_c_all,2)));
            if day_c_rel
                if isnan(C.data(day_c_all(day_c_rel)-1,1))
                    Y.Region(Reg).subregion(s).cases(i) =0;
                else
                    Y.Region(Reg).subregion(s).cases(i)   = C.data(day_c_all(day_c_rel(1))-1,1);% 1 is the UTLA
                end
            end
        end

        code_ages = find(agesperLA.Code == Y.Region(Reg).subregion(s).code);

        if( Reg ==6 && s == 2)
                code_ages = find(agesperLA.Code == 6000028) ;% Bournemouth:600028 % Christchurch:07000048 % Poole:06000029
               
 
        end
        
         if ( Reg ==6 && s == 6)
                code_ages = find(agesperLA.Code == 10000009) ; % Dorset
         end
     
        
        code_ages = code_ages(1);
        %----- ages in LA
        if ( Reg ==6 && s == 2)
             Y.Region(Reg).subregion(s).ages_0_90 =   ages_Bournemouth_Christhchurh_Poole(2:end);
        else
        for ages = 1:91
            Y.Region(Reg).subregion(s).ages_0_90(ages) =  eval(['agesperLA.VarName',num2str(ages+4),'(',num2str(code_ages),')']);
        end 
        end

        %---- dwellings (Rin) in regions - pop density (Rou) in regions
        reg_dwel_name = region_dwelling_table{Reg};
        
        dwelling_density =  eval(['dwellingnumberspopdensity.',reg_dwel_name]);
        all_dwelling     =  dwelling_density([1   4   7   8   9   10   12  13  14]); % one person, two or more unrelated adults, couple, couple 1-2 children ...see table
        tot_dwelling     =  sum( dwelling_density([1   4   7   8   9   10   12  13  14]));
        weightd          =  all_dwelling./tot_dwelling;
        Rin              =  [1  3  2  4  5   4    3    3   6]*weightd;
        Rou_density      =  dwelling_density(end); % people per square hectare
       
        Y.Region(Reg).subregion(s).Rin         = Rin;
        Y.Region(Reg).Rin                      = Rin;
        Y.Region(Reg).subregion(s).Rou_density = Rou_density;
        Y.Region(Reg).Rou_density              = Rou_density;
    end
    
    Y.Region(Reg).ages_0_90  = zeros(91,1);
    for  ages = 1:91
        for s = 1:length(subregi)
            Y.Region(Reg).ages_0_90(ages) = Y.Region(Reg).ages_0_90(ages) + Y.Region(Reg).subregion(s).ages_0_90(ages);
        end
    end
end


  for i = 1:9
        %---- dwellings (Rin) in regions - pop density (Rou) in regions
        reg_dwel_name = region_dwelling_table{i};
        
        dwelling_density =  eval(['dwellingnumberspopdensity.',reg_dwel_name]);
        all_dwelling     =  dwelling_density([1   4   7   8   9   10   12  13  14]); % one person, two or more unrelated adults, couple, couple 1-2 children ...see table
        tot_dwelling     =  sum( dwelling_density([1   4   7   8   9   10   12  13  14]));
        weightd          =  all_dwelling./tot_dwelling;
        Rin              =   [1  3  2   4   5   4    3    3   6]*weightd;
        Rin_for_NAT(i)              =   Rin ;
        Rou_density_for_NAT(i)      =  dwelling_density(end); % people per square hectare
        popR(i)                     =  Y.Region(i).pop;
  end
  
   Y.Nation.Rin                = Rin_for_NAT*(popR/(sum(popR)))';
   Y.Nation.Rou_density        = Rou_density_for_NAT*(popR/(sum(popR)))';
  
   
Y.dates = all_dates;
Y_prelim = Y;

region_ordered_nice = {'London', 'Yorkshire and The Humber','East Midlands',  'West Midlands' , 'East',...
    'South West', 'South East', 'North West', 'North East'}      ;

if plot_prelims
    
    %------- Plot National Statistics
    figure
    subplot(1,2,1)
    plot(Y_prelim.Nation.cases)
    tmp          = cumsum(Y_prelim.Nation.cases);
    total_cases  = tmp(end);
    title(['Total Cases as of ',all_dates{end},': ', num2str(total_cases)])
    
    subplot(1,2,2)
    plot(cumsum(Y_prelim.Nation.deaths))
    tmp          = cumsum(Y_prelim.Nation.deaths);
    total_deaths = tmp(end);
    title(['Total Deaths as of ',all_dates{end},' : ', num2str(total_deaths)])
     
    %------- Plot Regional Statistics
    for i = 1:9
        regional_all_cases(i,:) =   Y_prelim.Region(i).cases;
    end
    all_reg = sum(regional_all_cases,1);


    figure
    subplot(1,2,1)
    plot(Y_prelim.Nation.cases)
    hold on
    plot(all_reg ,'--r')
    subplot(1,2,2)
    plot(regional_all_cases')
    legend(region_ordered_nice) 
    
    %------- Plot Local Authority Statistics
    for i = 1:9
       
        regional_all_cases(i,:) =   Y_prelim.Region(i).cases;
        tmp = cumsum(Y_prelim.Region(i).cases);
        total_cases_region = tmp(end);
         LA_legend = cell(1,length(Y_prelim.Region(i).subregion));
         for j = 1:length(Y_prelim.Region(i).subregion)
             subregion_cases(j,:)   =  Y_prelim.Region(i).subregion(j).cases;
             tmp                    = cumsum(Y_prelim.Region(i).subregion(j).cases);
             tot_cases_subR(j)      = tmp(end);
             LA_legend{j}           = Y.Region(i).subregion(j).name{1};
         end
         
         all_subregion = sum(subregion_cases,1);
         
         figure
         subplot(1,2,1)
         plot(regional_all_cases(i,:))
         hold on
         plot( all_subregion,'--r')
         title(['Total Cases :' , region_ordered_nice{i}, ': ',  num2str(total_cases_region)])
         subplot(1,2,2)
         [x y] = sort(tot_cases_subR,'descend');
         plot(subregion_cases(y,:)')
         legend( LA_legend(y(1:8))  )
         
         for j = 1:8
          Y_prelim.Region(i).subregion_sorted(j).cases =  subregion_cases(y(j),:);
          Y_prelim.Region(i).subregion_sorted(j).LA    =  LA_legend(y(j));
          
         end
         
         clear LA_legend tot_cases_subR subregion_cases tmp
    end
  

end