(* ::Package:: *)

%%  This section computes VARs with different lag-lengths, and assesses the log-likelihoods, AIC, and BIC

clear; 
clc; 
for eyear=[1988 2015]
clear data

% Import Data and transform GNP into logs & difference 
gnp = xlsread ('gnp.xls'); 
urate = xlsread ('unrate.xls'); 

gnpl = log(gnp(:,1));
growth = 100*diff(gnpl(:,1)); 
data(:,1) = growth; 
data(:,2) = urate(2:end);

% Reporting dimensions

n = size(data,2);
syear = 1950;
data = data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

% Detrend unemployment

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:length(data),1)); 
data(1:95,1)= data(1:95,1)-mean1;
data(96:length(data),1)= data(96:length(data),1)-mean2;

% Asses information criteria for different VAR lenghts

for L= 2:1:16;
Mdl = varm(2,L);
EstMdl = estimate(Mdl,[data(:,1) data(:,2)]);
summarize(EstMdl)

end
end
