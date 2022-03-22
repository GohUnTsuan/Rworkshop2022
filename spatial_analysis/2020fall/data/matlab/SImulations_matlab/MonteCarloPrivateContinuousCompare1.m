% This file is made for the Monte Carlo experiment of the model of social
% interactions with private information. Consider the case that some
% characteristics are known only to an agent herself and that
% characteristics are discretely distributed.  

% The quadrature method is used to derive (numerically) the equilibrium
% conditional expectation functions.

% Codes are updated by Chao Yang, June 2019.

% REFERENCE:
% Yang and Lee, 2017, "Social interactions under incomplete information with
% heterogeneous expectations."?Journal of Econometrics?198.1 (2017): 65-83.

% First, input the simulated sample    
clear                            
tic;
initime          = cputime;             
time1            = clock;
rng(20130622);
load('MonteLinearPrivateContinuousDGP1_1.mat')                                                                    
[quadx,quadw]   = lgwt(K,-1,1);                                            % derive quadrature points and weights
global ExpM0 ExpMMis0
ExpM0           = zeros(K*n*G,1);
ExpMMis0        = zeros(n*G,1);                                          % set initial value of expectation
iter1           = 50;
TolX1           = 1.0000e-6; 
iter2           = 200;
TolX2           = 1.0000e-6; 
iter3           = 2500;                   
TolF            = 1.0000e-5;
TolX3           = 1.0000e-5;                                                                                                      

% Case 1 continuous model with public information   
% theta=[beta';log((1+lambda)/(1-lambda));sigma]';
% sigma=abs(theta(5));
% sigma^2=(theta(5))^2;
% Since |lambda|<1, we make the transformation that 
% lambda=(exp(theta(4))-1)/(exp(theta(4))+1).
% So theta(4)=log((1+lambda)/(1-lambda)). 
% Estimate mu, eta and rho in the first step using data about privately
% known characteristics (observed by econometricians)
% theta2=[mu;log((1+rho)/(1-rho))]
% theta2 is estimated using MLE
% Estimate of eta^2 is derived by first order condition

thetaXpest       = zeros(3,L);

recordthetaini_1 = zeros(5,L);                                             % store the initial guess by regression without network effects
thetaini_2       = 0.5*rand(5,1);                                          % initial guess randomly chosen
thetaini_3       = [0 0 0 0 1]';                                           % initial guess using 0 and 1  

thetaEst         = repmat([0 0 0 0.5 0.8]',1,L);   
logLEst          = ones(L,1);
exitcond         = 8*ones(L,1);
select           = zeros(L,1);
   
countloop        = 0; 
for l=1:L    
    W            = reshape(subW(:,l),n,n*G)';
    UseDDW       = sparse(reshape(repmat(1:(n*G),n,1),1,n*n*G),reshape(repmat(reshape(1:(n*G),n,G),n,1),1,n*n*G),...
         reshape(W',n*n*G,1),n*G,n*G);
    UseDDW       = UseDDW/friendnum; 
    UseDW        = sparse(reshape(repmat(1:(n*K*G),n,1),1,n*n*K*G),reshape(repmat(reshape(1:(n*K*G),n,K*G),n,1),1,n*n*K*G),...
        reshape(repmat(reshape(W',n*n,G),K,1),1,n*n*K*G),n*K*G,n*K*G);
    UseDW        = UseDW/friendnum;
    UseRDW       = sparse(reshape(repmat(1:(n*K*G),n,1),1,n*n*K*G),reshape(repmat(reshape(1:(n*K*G),n*K,G),n,1),1,n*n*K*G),...  
        reshape(repmat(reshape(W',n,n*G),K,1),1,n*K*n*G),n*K*G,n*K*G);
    UseRDW       = UseRDW/friendnum;
    UseXc        = Xc(:,l);
    UseXp        = Xp(:,l);                                                         % select sample used in the l-th repetition
    Usey_1_true  = y_1_true(:,l);
    
    n0            = n;
    G0            = G;
    FXPiter0      = iter1;
    Tol0          = TolX1;
    options       = optimset('MaxIter', iter3, 'TolFun', TolF, 'TolX', TolX3);  
    
    % estimate parameters for the distribution of Xp  
    theta2ini     = [0.5;1;1];
    theta2est     = fminunc(@(theta2) xplogL(theta2,UseXp,n0,G0),theta2ini,options);  
    muest         = theta2est(1);
    etasqest      = theta2est(2)^2;
    rhoest        = (exp(theta2est(3))-1)/(exp(theta2est(3))+1);  
    thetaXpest(:,l) = [muest;etasqest;rhoest];
    
    % Estimate the model without social interaction to derive one initial
    % guess
    UseX         = [ones(1,n*G);UseXc';UseXp'];
    if abs(det(UseX*(UseX')))>1.0000e-12
        betaini  = (UseX*(UseX'))\(UseX*Usey_1_true);
        sigmaini = sqrt(((Usey_1_true-(UseX')*betaini)')*(Usey_1_true-(UseX')*betaini)/(n*G));
        thetaini_1 = [betaini;0;sigmaini];  
    else
        thetaini_1 = [0.1 0.9 1.2 0 1.3]';
    end
    
    recordthetaini_1(:,l) = thetaini_1;  
    
    % estimation under true information structure
    mu               = muest;
    etasq            = etasqest;
    rho              = rhoest;
    K0               = K;
    
    theta2select=zeros(5,3);
    logL2select=zeros(3,1); 
    exit2select=8*ones(3,1); 
    [theta2select(:,1),logL2select(1),exit2select(1)]=fminsearch(@(theta) continuousPrivatecondlogL1...
        (theta,mu,etasq,rho,UseDW,UseRDW,UseXc,UseXp,Usey_1_true,n0,G0,K0,quadx,quadw,FXPiter0,Tol0),thetaini_1,options);
    [theta2select(:,2),logL2select(2),exit2select(2)]=fminsearch(@(theta) continuousPrivatecondlogL1...
        (theta,mu,etasq,rho,UseDW,UseRDW,UseXc,UseXp,Usey_1_true,n0,G0,K0,quadx,quadw,FXPiter0,Tol0),thetaini_2,options);
    [theta2select(:,3),logL2select(3),exit2select(3)]=fminsearch(@(theta) continuousPrivatecondlogL1...
        (theta,mu,etasq,rho,UseDW,UseRDW,UseXc,UseXp,Usey_1_true,n0,G0,K0,quadx,quadw,FXPiter0,Tol0),thetaini_3,options);
    [sortlogL,ranklogL]=sort(logL2select,'ascend');  
    thetaEst(:,l)=theta2select(:,ranklogL(1)); 
    logLEst(l)=-sortlogL(1);   
    exitcond(l)=exit2select(ranklogL(1));  
    select(l)=ranklogL(1);    
    
    
    countloop=countloop+1;                                                 % count the number of loops calculated
    
    save MonteLinearPrivateContinnuousCompareSmallK.mat L G n K Xc Xp subW y_1_true beta_true lambda_true sigmasq_true mu_true etasq_true rho_true...
    thetaXpest recordthetaini_1 thetaini_2 thetaini_3 thetaEst logLEst exitcond select countloop      
end

fintime = cputime;
elapsed = toc;
time2   = clock;
fprintf('TIC TOC: %g\n', elapsed);
fprintf('CPUTIME: %g\n', fintime - initime);
fprintf('CLOCK:   %g\n', etime(time2, time1));   
Ctime=fintime-initime; 
save MonteLinearPrivateContinuous1.mat L G n K Xc Xp subW y_1_true...
    beta_true lambda_true sigmasq_true mu_true etasq_true rho_true...
    thetaXpest recordthetaini_1 thetaini_2 thetaini_3...
    thetaEst logLEst exitcond select countloop... 
    elapsed Ctime friendnum iter1 iter2 TolX1 TolX2 iter3 TolF TolX3                                                                                                                                                                                                                                                                                                       




  

















