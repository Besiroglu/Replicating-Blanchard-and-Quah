(* ::Package:: *)

% This section replicates BQ's main results

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

% Demean and detrend unemployment

data(:,2)= data(:,2)-ones(length(data(:,2)),1)*mean(data(:,2)); 

data(:,2)=detrend(data(1:end,2));

% impose break in output growth

mean1=mean(data(1:95,1));
mean2=mean(data(96:151,1)); 

data(1:95,1)= data(1:95,1)-mean1;
data(96:151,1)= data(96:151,1)-mean2;

% Construct lags

L = 8;
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

% Configuring shocks
 
shockss = zeros (n,1); 
shocksd = zeros (n,1);

shockss(1) = 1; 
shocksd(2) = -1;

xs = zeros(n*L,1);
xs(1:n,1) = Q*shockss;

xd = zeros(n*L,1);
xd(1:n,1) = Q*shocksd;

% Setting up lag-matrices

T=40;

for t = 2:T
    xs(:,t) = F*xs(:,t-1); 
end 

for t = 2:T
    xd(:,t) = F*xd(:,t-1); 
end

% Defining paths response paths
% We employ cumulative sum as we're interested in LR effects

gnps = cumsum(xs (1,:)');
urates = xs (2,:)'; 

gnpd = cumsum(xd (1,:)'); 
urated = xd (2,:)';

% Bootstrapping at 95% with 1,000 replications

for i = 1:1000

    eb = e(ceil(nobs*rand(nobs,1)),:);
    
    % Initialize new data 
    
    Xn(1,:) = X(1,:); 
    Yn(1,:) = Xn(1,:)*beta+eb(1,:);
    
for  t=2:nobs
     Xn(t,:) = [Yn(t-1,:) Xn(t-1,1:end-n)]; 
     Yn(t,:) = Xn(t,:)*beta+eb(t,:); 
 end 
    
% Rerun estimation for each new pseudohistory 

betan = Xn\Yn; 
en = Yn - Xn*betan; 
Sigman = cov(en); 
F = companion (betan(1:end,:)); 

% Find new Q

I=eye(L*n,L*n);
I1=inv((I-F));
L1=[0 0; 0 0];
L1(1:2,1)=I1(1:2,1);
L1(1,1:2)=I1(1,1:2);
L1(2,2)=I1(2,2);
D=chol((L1*Sigma*L1'),'lower');
Q=inv(L1)*D;

% Initialise shocks \[Dash] as usual

xs = zeros(n*L,1); 
xd = zeros(n*L,1); 

xs(1:n,1) = Q*shockss;
xd(1:n,1) = Q*shocksd; 
    
   for t = 2:T
    xs(:,t) = F*xs(:,t-1); 
    xd(:,t) = F*xd(:,t-1);
       
end

   
   gnpsn(:,i)=cumsum(xs (1,:)');
   uratesn(:,i)=xs (2,:)'; 
   
   gnpdn(:,i)=cumsum(xd (1,:)');
   uratedn(:,i)=xd (2,:)'; 
end

	% prepare IRFs for plotting 

conf = 0.95; 
xaxis = (1:T)'; 

% Supply shock confidence bounds 
gnphs = quantile (gnpsn',conf+(1-conf)/2)'; 
gnpls = quantile (gnpsn',(1-conf)/2)'; 
gnpms = median(gnpsn,2); 

uratehs = quantile (uratesn',conf+(1-conf)/2)';
uratels = quantile (uratesn',(1-conf)/2)';
uratems = median(uratesn,2); 

% Demand shock confidence bounds 
gnphd = quantile (gnpdn',conf+(1-conf)/2)'; 
gnpld = quantile (gnpdn',(1-conf)/2)'; 
gnpmd = median(gnpdn,2); 

uratehd = quantile (uratedn',conf+(1-conf)/2)';
urateld = quantile (uratedn',(1-conf)/2)';
uratemd = median(uratedn,2); 

% Plot

figure (1); 

subplot(2,2,1)
[hl,hp] = boundedline(xaxis,gnpms,[gnpms-gnpls gnphs-gnpms],'alpha'); 
box on 
set (hl,'linewidth',1.6) 
title('Supply shock - GNP growth','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex') 

subplot(2,2,2)
[hl,hp] = boundedline(xaxis,uratems,[uratems-uratels uratehs-uratems],'alpha'); 
box on 
set (hl,'linewidth',1.6) 
title('Supply shock - unemployment rate','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex') 

subplot(2,2,3)
[hl,hp] = boundedline(xaxis,gnpmd,[gnpmd-gnpld gnphd-gnpmd],'alpha'); 
box on 
set (hl,'linewidth',1.6) 
title('Demand shock - GNP growth','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex') 

subplot(2,2,4)
[hl,hp] = boundedline(xaxis,uratemd,[uratemd-urateld uratehd-uratemd],'alpha'); 
box on 
set (hl,'linewidth',1.6) 
title('Demand shock - unemployment rate','Fontsize',16,'fontname','arial','Interpreter','latex') 
set(gca,'Fontsize', 14, 'fontname','arial')
ylabel ('Deviation (\%)','Fontsize', 14, 'fontname','arial','Interpreter','latex')
xlabel ('Time (quarters)','Fontsize',14,'fontname','arial','Interpreter','latex')
