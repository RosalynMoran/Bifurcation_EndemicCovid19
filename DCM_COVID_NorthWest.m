% adapted for England, Regions and Local Authorities from
% FORMAT [DCM] = DEM_COVID_X(data)
% data    - data    to model [default: data = DATA_COVID_JHU]
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine illustrates the Bayesian model inversion of a generative
% model of coronavirus spread using variational techniques (variational
% Laplace). This (pandemic) model is composed of regional (epidemic)
% models. In brief, the model for a single region comprises four factors,
% each with four states, giving 256 states or compartments per region.
% These regional models are then assembled to model the coupling among
% eight regions giving 256^8 compartments. However, due to certain
% conditional independencies, this can be treated as a collection of 256
% compartmental models; providing one carefully links the state of one
% region to the state of another. Here, this linking or connectivity is
% parameterised in terms of a probability flux or exchange of people from
% one regional population to another. Regional factors include location,
% immune status, clinical status and testing status. The transitions among
% states of any factor depends upon other factors. For example, the
% probability that I will move from a state of being asymptomatic to being
% symptomatic depends upon whether I am infected or not. Similarly, the
% probability that I will move from one region to another depends upon
% whether I am at work (i.e., not at home). In short, the exchange between
% different regional populations is limited to the people who are not at
% home and are consequently in a position to travel. The parameters of
% interregional coupling correspond to rate constants or effective
% connectivity that can be reciprocal and asymmetric. For example, the
% probability of moving to New York from New Jersey does not have to be the
% same as a probability of moving from New Jersey to New York. Note that
% the movement between regions can be restricted to a chain. In other
% words, to get from the first state to the last state, I have to go
% through all other states.
%
% Each subsection produces one or two figures that are described in the
% annotated (Matlab) code. These subsections call various subroutines that
% provide a more detailed description of things like the generative model,
% its priors and the evaluation of confidence intervals.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_X.m 7838 2020-04-23 17:40:45Z karl $
clear all
close all
 
% set path spm folder
addpath(genpath('C:\Users\rosal\Dropbox\Covid_19\spm12_latest\spm12'))
 

% Invert England ----  Deaths - Cases
%==========================================================================


% -----------------------------------1.   Get Data

plot_prelims =   0;   % plot to return subregions sorted
Y_prelim     =   DATA_COVID_UK_Prelims(plot_prelims);
Nation       =   0;
Region       =   'all'; % or insert numbers 1 -9 for regional local authority network
regioni      =   8;   % The Northwest
[Y , pop]    =   DATA_COVID_UK_PHE(Y_prelim, Nation, Region);



% -----------------------------------2.  Set Up Priors

[pE,pC,str] = spm_COVID_priors;        % priors
pE.N        = log(pop(regioni));                % fix (log) population size in Millions  
pC.N        = 0;                       % with no uncertainty
pE.r        = log(1/4);                % assume 25% resistance
pC.r        = 1/64;                    % with some uncertainty



%-----------------------------------3.  Do inversion based on data Y and priors: Y(:,1) = Daily Death Reports, Y(:,2) = Daily Case Reports

[F,Ep_Region,Cp_Region] = spm_COVID_cases_only(Y(:,:,regioni),pE,pC);

DCM.F = F;
DCM.Ep = Ep_Region;
DCM.Cp = Cp_Region;
DCM.Data   = Y;

%-----------------------------------4.  Get Model Fits or Forecast (using
% M.T to extend days into future,e.g. Y is length 99 days, M.T = 199 to
% predict 100 days into the future

M.T                        =   length(Y);
U                          =   [1:6];     % all outputs Y(:,1):deaths, Y(:,2):cases, Y(:,3):CCUOccupancy, Y(:,4):EffectiveReproductionRate(R)  Y(:,5):immunity  Y(:,6):total number of daily tests (positive or negative)

[YRegionPred ,XRegionPred] =   spm_COVID_gen(Ep_Region,M,1:6);

DCM.R0     = YRegionPred(:,4);
DCM.Latent_States = XRegionPred;

%-----------------------------------5.   Plot things

figure
subplot(2,2,1)
plot(XRegionPred{2})
legend({'Susceptible','Infected','Infectious','Immune','Resistant'})
 
subplot(2,2,2)
plot(YRegionPred(:,2))
hold on
plot(Y(:,:,regioni))
title('Daily Reported Cases')

subplot(2,2,3)
plot(cumsum(YRegionPred(:,1)))
title('Deaths, cumulative')
 
subplot(2,2,4)
plot(YRegionPred(:,4),'k')
title('R_t')



%--------------------------------------------------------6.  simulate 40% FTTI with 16 month immunity

% ----- Initial Conditions

s        =  XRegionPred{2}(length(Y),1); % susceptible
n        =  XRegionPred{2}(length(Y),3); % initial infected
m        =  XRegionPred{2}(length(Y),4); % immune
r        =  XRegionPred{2}(length(Y),5); % resistant
s = s+r-0.2;
r = 0.2;

% ----- Key Parameters

Ep_Region.Tim                     =  log(16);      % months of immunity
Ep_Region.ttt                     =  log(0.90);    % ttt 10% effective
M.T                               =  30*24;        % days from may 20th to simulate


[YPredFut ,XPredFut]  =  spm_COVID_gen_ICs(Ep_Region,M,U,n,m,r,s);


%-----------------------------------7.   Plot things in the futute

figure
subplot(2,2,1)
plot(YPredFut(:,2))
title('Daily Reported Cases, Future predictions')

subplot(2,2,2)
plot(cumsum(YPredFut(:,1)))
title('Deaths, Cumulative, Future predictions')
 
subplot(2,2,3)
plot( (YPredFut(:,1)))
title('Daily Deaths , Future predictions')

subplot(2,2,4)
plot(YPredFut(:,6) ) 
title('Tests needed, Future predictions') 
 
 


%   Inversion for Single Time Series of Positive Tests from Regions and Local
%   Authorities
% 
%   ----------------------------------- Invert Regional Network  -----  Cases Only
%
% Nation       =    0;
% Region       =    'all'; % or inserts numbers 1 -9 for regional local authority network
% [YR , popR]  =    DATA_COVID_UK_PHE(Y_prelim, Nation, Region);
% pE           =    Ep_Nation;
% pC           =    spm_unvec(diag(Cp_Nation),pC);
% 
% for i = 1:9
%     
%  
%     pE.N       =    log(popR(i));                % fix (log) population size
%     pC.N       =    0;                       % with no uncertainty
%     pE.r       =    log(1/4);                % assume 25% resistance
%     pC.r       =    1/64;
%     
%     [F,Ep,Cp] = spm_COVID_cases_only(YR(:,:,i),pE,pC);
%     
%     % assemble prior and posterior estimates (and log evidence)
%     %------------------------------------------------------------------
%     M.Q(i).Ep  = Ep;
%     M.Q(i).Cp  = Cp;
%     M.Q(i).F   = F;
%     M.Q(i).YR  = YR;
%     
% end
% 
 
 

