(* ::Package:: *)

%%  Calculating optimal lags in the VAR model for different time periods

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

%% Lagrange Multiplier test Ljung-Box Q-test for different models

clear; 
clc; 

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
eyear = 1988;

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

% Construct lags

for L =[1 2 4 6 8 10 12 14 16]

LM = Lmatrix(data,0:L);
n=size(data,2);
X= LM(L+1:end,n+1:end);
Y = LM(L+1:end,1:n);
beta=X\Y;

e= Y-X*beta;

% Lagrange Multiplier (LM) test for autoregressive conditional heteroscedasticity (archtest)

archtest(e(:,1),'Alpha',0.05)
archtest(e(:,2),'Alpha',0.05)

% Ljung-Box Q-test for residual autocorrelation (lbqtest)

lbqtest(e(:,1),'Alpha',0.05)
lbqtest(e(:,2),'Alpha',0.05)

end
