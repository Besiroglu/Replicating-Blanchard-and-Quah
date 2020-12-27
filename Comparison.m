(* ::Package:: *)

%% Evaluating the effect of supply disturbances on unemployment

clear data
clear ys

% Import Data and transform GNP into logs & difference 
gnp = xlsread ('gnp.xls'); 
urate = xlsread ('unrate.xls'); 

gnpl = log(gnp(:,1));

growth = 100*diff(gnpl(:,1)); 

data(:,1) = growth; 
data(:,2) = urate(2:end);

% Reporting dimensions

n = size(data,2);

eyear=2015;
syear = 1950;

data = data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:length(data),1)); 

data(1:95,1)= data(1:95,1)-mean1;
data(96:length(data),1)= data(96:length(data),1)-mean2;

% Vary the number of lags between 2 and 8

for L= 2:2:8;

LM = Lmatrix(data,0:L);
n=size(data,2);
X= LM(L+1:end,n+1:end);
Y = LM(L+1:end,1:n);
beta=X\Y;

nobs=size(Y,1);

% define companion matrix

F= companion(beta(1:end,:));

% calculate Sigma to find D,Q
e= Y-X*beta;
Sigma = cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;
% sanity check: Q*Q' = Sigma

% Set shocks  (supply shocks)

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
ys(:,i) = cumsum(x (1,:)');
uns(:,i)= x (2,:)';
end
lowestunemployment1(:,L/2)=(max(uns(:,:)))';
end

eyear=1988;

data = data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:length(data),1)); 

data(1:95,1)= data(1:95,1)-mean1;
data(96:length(data),1)= data(96:length(data),1)-mean2;

% Vary the number of lags between 2 and 8

for L= 2:2:8;

LM = Lmatrix(data,0:L);
n=size(data,2);
X= LM(L+1:end,n+1:end);
Y = LM(L+1:end,1:n);
beta=X\Y;

nobs=size(Y,1);

% define companion matrix

F= companion(beta(1:end,:));

% calculate Sigma to find D,Q
e= Y-X*beta;
Sigma = cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;
% sanity check: Q*Q' = Sigma

% Set shocks  (supply shocks)

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
ys(:,i) = cumsum(x (1,:)');
uns(:,i)= x (2,:)';
end
lowestunemployment2(:,L/2)=(max(uns(:,:)))';
end

figure (4);
subplot(1,2,1)
boxplot(lowestunemployment1)
ylim([-0.4 0.8])
title('Highest level of unemployment due to supply shock (1950-2015)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex')
set(gca,'XTickLabel',[2:2:8] );

subplot(1,2,2)
boxplot(lowestunemployment2)
ylim([-0.4 0.8])
title('Highest level of unemployment due to supply shock (1950-1988)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex')
set(gca,'XTickLabel',[2:2:8] );


%% Evaluating the long run effect of supply disturbances

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

eyear=1988;
syear = 1950;

data = data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:length(data),1)); 

data(1:95,1)= data(1:95,1)-mean1;
data(96:length(data),1)= data(96:length(data),1)-mean2;

% Vary the number of lags between 2 and 8

for L= 2:2:8;

LM = Lmatrix(data,0:L);
n=size(data,2);
X= LM(L+1:end,n+1:end);
Y = LM(L+1:end,1:n);
beta=X\Y;

nobs=size(Y,1);

% define companion matrix

F= companion(beta(1:end,:));

% calculate Sigma to find D,Q
e= Y-X*beta;
Sigma = cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;
% sanity check: Q*Q' = Sigma

% Set shocks  (supply shocks)

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
ys(:,i) = cumsum(x (1,:)');
uns(:,i)= x (2,:)';
end

Longruneffect1(:,L/2)=ys(end,:);

end

%% VAR for period 1950-2015

clear data
clear ys

% Import Data and transform GNP into logs & difference 
gnp = xlsread ('gnp.xls'); 
urate = xlsread ('unrate.xls'); 

gnpl = log(gnp(:,1));

growth = 100*diff(gnpl(:,1)); 

data(:,1) = growth; 
data(:,2) = urate(2:end);

% Reporting dimensions

n = size(data,2);

eyear=2015;
syear = 1950;

data = data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:length(data),1)); 

data(1:95,1)= data(1:95,1)-mean1;
data(96:length(data),1)= data(96:length(data),1)-mean2;

% Vary the number of lags between 2 and 8

for L= 2:2:8;

LM = Lmatrix(data,0:L);
n=size(data,2);
X= LM(L+1:end,n+1:end);
Y = LM(L+1:end,1:n);
beta=X\Y;

nobs=size(Y,1);

% define companion matrix

F= companion(beta(1:end,:));

% calculate Sigma to find D,Q
e= Y-X*beta;
Sigma = cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;
% sanity check: Q*Q' = Sigma

% Set shocks  (supply shocks)

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
ys(:,i) = cumsum(x (1,:)');
uns(:,i)= x (2,:)';
end
Longruneffect2(:,L/2)=ys(end,:);
end

figure (3);
subplot(1,2,1)
boxplot(Longruneffect1)
ylim([-0.5 2.5])
title('Long run effect of supply shock on GNP (1950-1988)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex')
set(gca,'XTickLabel',[2:2:8] );

subplot(1,2,2)
boxplot(Longruneffect2)
ylim([-0.5 2.5])
title('Long run effect of supply shock on GNP (1950-2015)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex')
set(gca,'XTickLabel',[2:2:8] );


%% half lifes for different lags

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
eyear = 2015;

data = data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:151,1)); 

data(1:95,1)= data(1:95,1)-mean1;
data(96:151,1)= data(96:151,1)-mean2;

% Vary the number of lags between 2 and 8

for L= 2:2:8;

LM = Lmatrix(data,0:L);
n=size(data,2);
X= LM(L+1:end,n+1:end);
Y = LM(L+1:end,1:n);
beta=X\Y;

nobs=size(Y,1);

% define companion matrix

F= companion(beta(1:end,:));

% calculate Sigma to find D,Q
e= Y-X*beta;
Sigma = cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;
% sanity check: Q*Q' = Sigma

% Set shocks  (supply shocks)

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
ys(:,i) = cumsum(x (1,:)');
uns(:,i)= x (2,:)';
end

x=(1:T)';

% Counter 
Co=L/2;

% Set shocks  (demand shocks)

shocks = zeros(n,1);
shocks (2)= -1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (2)= -1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
yd(:,i) = cumsum(x (1,:)');
und(:,i)= x (2,:)';
end

% Calculate half-lifes for each response functions, store these in matrix

for j =1:1000;
 try
lower1=(find(und(1:T,j)==min(und(1:T,j)),j,'first'));
halflifeud(j,L/2)=find(und(lower1:T,j)>=0.5*min(und(1:T,j)),1,'first');
    % catch if error
  catch
    fprintf('Inconsistent data in iteration % s, skipped.\n', j);
  end
end 

for j =1:1000;
 try
lower2=find(yd(1:T,j)==max(yd(1:T,j)),j,'first');
halflifeyd(j,L/2)=find(yd(lower2:T,j)<=0.5*max(yd(1:T,j)),1,'first');
    % catch if error
  catch
    fprintf('Inconsistent data in iteration % s, skipped.\n', j);
  end
end 
end

eyear = 1988;

data = data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:151,1)); 

data(1:95,1)= data(1:95,1)-mean1;
data(96:151,1)= data(96:151,1)-mean2;

% Vary the number of lags between 2 and 8

for L= 2:2:8;

LM = Lmatrix(data,0:L);
n=size(data,2);
X= LM(L+1:end,n+1:end);
Y = LM(L+1:end,1:n);
beta=X\Y;

nobs=size(Y,1);

% define companion matrix

F= companion(beta(1:end,:));

% calculate Sigma to find D,Q
e= Y-X*beta;
Sigma = cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;
% sanity check: Q*Q' = Sigma

% Set shocks  (supply shocks)

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (1)= 1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
ys(:,i) = cumsum(x (1,:)');
uns(:,i)= x (2,:)';
end

x=(1:T)';

% Counter 
Co=L/2;

% Set shocks  (demand shocks)

shocks = zeros(n,1);
shocks (2)= -1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t = 2:T;
x(:,t)=F*x(:,t-1);
end

% Bootstrapping 1000 times

for i= 1:1000;

eb= e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;

betan = Xn\Yn;
en = Yn-Xn*betan;

Sigman= cov(en);
F= companion (betan (1:end,:));

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;

shocks = zeros(n,1);
shocks (2)= -1;
x= zeros(n*L,1);
x(1:n,1)=Q*shocks;

for t = 2:T;
x(:,t)=F*x(:,t-1);
end
yd(:,i) = cumsum(x (1,:)');
und(:,i)= x (2,:)';
end

% Calculate half-lifes for each response functions, store these in matrix

for j =1:1000;
 try
lower1=(find(und(1:T,j)==min(und(1:T,j)),j,'first'));
halflifeud2(j,L/2)=find(und(lower1:T,j)>=0.5*min(und(1:T,j)),1,'first');
    % catch if error
  catch
    fprintf('Inconsistent data in iteration % s, skipped.\n', j);
  end
end 

for j =1:1000;
 try
lower2=find(yd(1:T,j)==max(yd(1:T,j)),j,'first');
halflifeyd2(j,L/2)=find(yd(lower2:T,j)<=0.5*max(yd(1:T,j)),1,'first');
    % catch if error
  catch
    fprintf('Inconsistent data in iteration % s, skipped.\n', j);
  end
end 
end

%% create plots

figure (3);
subplot(2,2,1)
boxplot(halflifeyd)
ylim([1 15])
box on 
title('Half-life of demand shock on GNP (1950-2015)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Time (quarters)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex') 
set(gca,'XTickLabel',[2:2:8] );

subplot(2,2,2)
boxplot(halflifeud)
ylim([1 15])
title('Half-life of demand shock on unemployment (1950-2015)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Time (quarters)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex')
set(gca,'XTickLabel',[2:2:8] );

figure (3);
subplot(2,2,3)
boxplot(halflifeyd2)
ylim([1 15])
box on 
title('Half-life of demand shock on GNP (1950-1988)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Time (quarters)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex') 
set(gca,'XTickLabel',[2:2:8] );

subplot(2,2,4)
boxplot(halflifeud2)
ylim([1 15])
title('Half-life of demand shock on unemployment (1950-1988)','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Time (quarters)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Number of lags in VAR','Fontsize',14,'fontname','arial','Interpreter','latex')
set(gca,'XTickLabel',[2:2:8] );

