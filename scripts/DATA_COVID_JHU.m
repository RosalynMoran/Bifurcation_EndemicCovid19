function [Y] = DATA_COVID_JHU(country)

% Data retrieval function for COVID modelling
% Using JHU data
% CCSE
%
% https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
 

 
DD = importdata('time_series_covid19_deaths_global.csv');
CC = importdata('time_series_covid19_confirmed_global.csv');
 
country_index = find(strcmp(CC.textdata(:,2),'Iran'))-1;

cases  = CC.data(country_index,3:end);
deaths = DD.data(country_index,3:end);

Y(:,1)    = diff(deaths); % per day
Y(:,2)    = diff(cases);  % per day
