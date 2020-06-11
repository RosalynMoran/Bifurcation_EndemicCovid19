function [Y,pop] = DATA_COVID_UK_PHE(Y_prelim,Nation,Region)

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
 
% Y returns as [Nregion/subregions x 2 or 1 times Ndays]

% See Yin for data structure
%==========================================================================
% s     = 2;                                          % data smoothing (days)
% T     = 2;                                         % days to skip from {'2020-01-30'};
% k     = 1;
%  ----  using simple 3 day moving average  over daily deaths and cases
%  (England) or just cases (regionals & Local Authority)

% England  ----------------------------------- deaths - cases
if Nation == 1
    
    Y(:,1) =  movmean(Y_prelim.Nation.deaths,3) ;
    Y(:,2)  =  movmean(Y_prelim.Nation.cases,3) ;
 
%     YJHUD = gradient(spm_conv([zeros(8,1); Y_prelim.Nation.deaths],2))
%     YJHUC = gradient(spm_conv([zeros(8,1); Y_prelim.Nation.cases],2))
    
    % population  
    %----------------------------------------------------------------------
    pop = Y_prelim.Nation.pop;
end
    


% Regions -----------------------------------------  cases only
if Nation == 0
    if strcmp(Region,'all') 
        
        for i = 1:9
            
             Y(:,1,i) = movmean(Y_prelim.Region(i).cases,3);
             pop(i) = Y_prelim.Region(i).pop;
             
        end
    else
        
        no_subregions = length(Y_prelim.Region(Region).subregion);
        
        for i = 1:no_subregions
            
            Y(:,1,i) =  movmean(Y_prelim.Region(Region).subregion(i).cases,3);
            pop(i)   =  Y_prelim.Region(Region).subregion(i).pop*1e-6; % in millions
            
        end
    end
end

  

   
    

    

       
      
% %     % data
% %     %----------------------------------------------------------------------
% %     Y(:,1,i)    = Data(i).death;
% %     Y(:,2,i)    = Data(i).cases;
% %     
% % end
% % 
% % % remove negative values
% % %--------------------------------------------------------------------------
% % Y(Y < 0) = 0;
% % 
% % 
% % 
% % 
% 
