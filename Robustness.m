(* ::Package:: *)

%% This section checks the robustness of B&Q's results to the number of lags in the VAR models

clear;
clc;

% Import Data and transform GNP into logs & difference

gnp=xlsread('gnp.xls');
urate=xlsread('unrate.xls');
gnpl=log(gnp(:,1));
growth=100*diff(gnpl(:,1));
data(:,1)=growth;
data(:,2)=urate(2:end);

% Reporting dimensions

n=size(data,2);
syear=1950;
eyear=1988;
data=data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)=data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2));
data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:151,1));
data(1:95,1)=data(1:95,1)-mean1;
data(96:151,1)=data(96:151,1)-mean2;

% Vary the number of lags between 4 and 16

for L=2:2:8;
LM=Lmatrix(data,0:L);
n=size(data,2);
X=LM(L+1:end,n+1:end);
Y=LM(L+1:end,1:n);
beta=X\Y;
nobs=size(Y,1);

% define companion matrix

F=companion(beta(1:end,:));

% calculate Sigma to find D,Q

e=Y-X*beta;
Sigma=cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[00;00];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;

% sanity check:Q*Q'=Sigma
% Set shocks(supply shocks)

shocks=zeros(n,1);
shocks(1)=1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
yp=cumsum(x (1,:)');
up=x (2,:)';

% Bootstrapping 1000 times

for i=1:1000;
eb=e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;
betan=Xn\Yn;
en=Yn-Xn*betan;
Sigman=cov(en);
F=companion(betan(1:end,:));
I=eye(L*n,L*n);
I1=inv((I-F));
L1=[00;00];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;
shocks=zeros(n,1);
shocks(1)=1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
y(:,i)=cumsum(x (1,:)');
un(:,i)=x (2,:)';
end
conf=0.95;
x=(1:T)';
yh=quantile (y',conf+(1-conf)/2)';
yl=quantile (y',(1-conf)/2)';
yms=median(y,2);
unh=quantile (un',conf+(1-conf)/2)';
unl=quantile (un',(1-conf)/2)';
unms=median(un,2);

% Create counter

Co=L/2;

% Store each median for supply-IRFs

Coys1(:,Co)=yms;
Cous1(:,Co)=unms;

% Set demand shocks

shocks=zeros(n,1);
shocks(2)=-1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
yp=cumsum(x (1,:)');
up=x (2,:)';

% Bootstrapping 1000 times

for i=1:1000;
eb=e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;
betan=Xn\Yn;
en=Yn-Xn*betan;
Sigman=cov(en);
F=companion(betan(1:end,:));
I=eye(L*n,L*n);
I1=inv((I-F));
L1=[00;00];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;
shocks=zeros(n,1);
shocks(2)=-1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
y(:,i)=cumsum(x (1,:)');
un(:,i)=x (2,:)';
end
conf=0.95;
x=(1:T)';
yh=quantile (y',conf+(1-conf)/2)';
yl=quantile (y',(1-conf)/2)';
ymd=median(y,2);
unh=quantile (un',conf+(1-conf)/2)';
unl=quantile (un',(1-conf)/2)';
unmd=median(un,2);
 
% Store each median for demand-IRFs

Coyd1(:,Co)=ymd;
Coud1(:,Co)=unmd;
end
clear data

% Import Data and transform GNP into logs & difference 

gnp=xlsread('gnp.xls');
urate=xlsread('unrate.xls');
gnpl=log(gnp(:,1));
growth=100*diff(gnpl(:,1));
data(:,1)=growth;
data(:,2)=urate(2:end);

% Reporting dimensions

n=size(data,2);
syear=1950;
eyear=2015;
data=data((syear-1950)*4+1:(eyear-1950)*4-1,:);

% Demean and detrend unemployment

data(:,2)=data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2));
data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:length(data),1));
data(1:95,1)=data(1:95,1)-mean1;
data(96:length(data),1)=data(96:length(data),1)-mean2;

% Vary the number of lags between 4 and 8

for L=2:2:8;
LM=Lmatrix(data,0:L);
n=size(data,2);
X=LM(L+1:end,n+1:end);
Y=LM(L+1:end,1:n);
beta=X\Y;
nobs=size(Y,1);

% define companion matrix

F=companion(beta(1:end,:));

% calculate Sigma to find D,Q

e=Y-X*beta;
Sigma=cov(e);

% Find Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[00;00];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;
% sanity \[LongDash] check that Q*Q'=Sigma

% Set shocks(supply shocks)

shocks=zeros(n,1);
shocks(1)=1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
yp=cumsum(x (1,:)');
up=x (2,:)';

% Bootstrapping 1000 times

for i=1:1000;
eb=e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;
betan=Xn\Yn;
en=Yn-Xn*betan;
Sigman=cov(en);
F=companion(betan(1:end,:));
I=eye(L*n,L*n);
I1=inv((I-F));
L1=[00;00];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;
shocks=zeros(n,1);
shocks(1)=1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
y(:,i)=cumsum(x (1,:)');
un(:,i)=x (2,:)';
end
conf=0.95;
x=(1:T)';
yh=quantile (y',conf+(1-conf)/2)';
yl=quantile (y',(1-conf)/2)';
yms=median(y,2);
unh=quantile (un',conf+(1-conf)/2)';
unl=quantile (un',(1-conf)/2)';
unms=median(un,2);

% Counter

Co=L/2;

% Store each median for supply-IRFs

Coys2(:,Co)=yms;
Cous2(:,Co)=unms;

% Set demand shocks

shocks=zeros(n,1);
shocks(2)=-1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;

% IRF replicating with 40 periods

T=40;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
yp=cumsum(x (1,:)');
up=x (2,:)';

% Bootstrapping 1000 times

for i=1:1000;
eb=e(ceil(nobs*rand(nobs,1)),:);
Xn=X;
Yn=Xn*beta+eb;
betan=Xn\Yn;
en=Yn-Xn*betan;
Sigman=cov(en);
F=companion(betan(1:end,:));
I=eye(L*n,L*n);
I1=inv((I-F));
L1=[00;00];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigman*L1'),'lower');
Q=inv(L1)*D;
shocks=zeros(n,1);
shocks(2)=-1;
x=zeros(n*L,1);
x(1:n,1)=Q*shocks;
for t=2:T;
x(:,t)=F*x(:,t-1);
end
y(:,i)=cumsum(x (1,:)');
un(:,i)=x (2,:)';
end
conf=0.95;
x=(1:T)';
yh=quantile (y',conf+(1-conf)/2)';
yl=quantile (y',(1-conf)/2)';
ymd=median(y,2);
unh=quantile (un',conf+(1-conf)/2)';
unl=quantile (un',(1-conf)/2)';
unmd=median(un,2);

% Store each median for demand-IRFs

Coyd2(:,Co)=ymd;
Coud2(:,Co)=unmd;
end
xaxis=(1:T)';
figure(2);
subplot(2,4,1)
plot(Coys1(:,1),'g','linewidth',1.5)
hold on;
plot(Coys1(:,2),'LineStyle','--','linewidth',1.5)
plot(Coys1(:,3),'LineStyle','--','linewidth',1.5)
plot(Coys1(:,4),'k','linewidth',1.5)
hold off
box on
title('Supply shock - GNP growth','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
subplot(2,4,2)
plot(Cous1(:,1),'g','linewidth',1.5)
hold on;
plot(Cous1(:,2),'LineStyle','--','linewidth',1.5)
plot(Cous1(:,3),'LineStyle','--','linewidth',1.5)
plot(Cous1(:,4),'k','linewidth',1.5)
hold off
box on
title('Supply shock - unemployment rate','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
subplot(2,4,3)
plot(Coys2(:,1),'g','linewidth',1.5)
hold on;
plot(Coys2(:,2),'LineStyle','--','linewidth',1.5)
plot(Coys2(:,3),'LineStyle','--','linewidth',1.5)
plot(Coys2(:,4),'k','linewidth',1.5)
hold off
box on
title('Supply shock - GNP growth','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
subplot(2,4,4)
plot(Cous2(:,1),'g','linewidth',1.5)
hold on;
plot(Cous2(:,2),'LineStyle','--','linewidth',1.5)
plot(Cous2(:,3),'LineStyle','--','linewidth',1.5)
plot(Cous2(:,4),'k','linewidth',1.5)
hold off
box on
title('Supply shock - unemployment rate','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
subplot(2,4,5)
plot(Coyd1(:,1),'g','linewidth',1.5)
hold on;
plot(Coyd1(:,2),'LineStyle','--','linewidth',1.5)
plot(Coyd1(:,3),'LineStyle','--','linewidth',1.5)
plot(Coyd1(:,4),'k','linewidth',1.5)
hold off
box on
title('Demand shock - GNP growth','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
subplot(2,4,6)
plot(Coud1(:,1),'g','linewidth',1.5)
hold on;
plot(Coud1(:,2),'LineStyle','--','linewidth',1.5)
plot(Coud1(:,3),'LineStyle','--','linewidth',1.5)
plot(Coud1(:,4),'k','linewidth',1.5)
hold off
box on
title('Demand shock - unemployment rate','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
subplot(2,4,7)
plot(Coyd2(:,1),'g','linewidth',1.5)
hold on;
plot(Coyd2(:,2),'LineStyle','--','linewidth',1.5)
plot(Coyd2(:,3),'LineStyle','--','linewidth',1.5)
plot(Coyd2(:,4),'k','linewidth',1.5)
hold off
box on
title('Demand shock - GNP growth','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
subplot(2,4,8)
plot(Coud2(:,1),'g','linewidth',1.5)
hold on;
plot(Coud2(:,2),'LineStyle','--','linewidth',1.5)
plot(Coud2(:,3),'LineStyle','--','linewidth',1.5)
plot(Coud2(:,4),'k','linewidth',1.5)
hold off
box on
title('Demand shock - unemployment rate','Fontsize',16,'fontname','arial','Interpreter','latex')
set(gca,'Fontsize',14,'fontname','arial')
ylabel('Deviation (\%)','Fontsize',14,'fontname','arial','Interpreter','latex')
xlabel('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
