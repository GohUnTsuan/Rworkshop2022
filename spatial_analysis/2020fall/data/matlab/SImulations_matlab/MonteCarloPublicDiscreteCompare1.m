% This file is made for the Monte Carlo experiment of the model of social
% interactions with private information. This set of experiments consider
% the case that all the variables available in the data sets are also
% publcly known to all group members.

% Codes are written and updated by Chao Yang.

% REFERENCE:
% Yang and Lee, 2017, "Social interactions under incomplete information with
% heterogeneous expectations."?Journal of Econometrics?198.1 (2017): 65-83


% First, input the simulated sample    
clear                            
tic;
initime          = cputime;            
time1            = clock;
rng(20130622); 
load('MonteContinuousPublicDiscreteCompare1.mat')                
global ExpM0 ExpMMis0
ExpM0           = zeros(n*G,1);
ExpMMis0        = zeros(2*n*G,1);                                          % set initial value of expectation
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

recordthetaini_1 = zeros(5,L);                                             % store the initial guess by regression without network effects
thetaini_2       = 0.5*rand(5,1);                                          % initial guess randomly chosen
thetaini_3       = [0 0 0 0 1]';                                           % initial guess using 0 and 1 
thetaini_4       = [0 0 0 log(1.99/0.01) 1]'; 

thetaEst         = repmat([0 0 0 0.5 0.8]',1,L);    
logLEst          = ones(L,1);
exitcond         = 8*ones(L,1);
select           = zeros(L,1);
    
countloop        = 0; 
for l=1:L    
    W            = Wrecord(:,l); 
    UseDDW       = sparse(reshape(repmat(1:(n*G),n,1),1,n*n*G),reshape(repmat(reshape(1:(n*G),n,G),n,1),1,n*n*G),...
         reshape(W',n*n*G,1),n*G,n*G);
    UseDDW       = UseDDW/friendnum; 
    UseXc        = Xc(:,l);
    UseXp        = Xp(:,l);                                                         % select sample used in the l-th repetition
    Usey_1_true  = yrecord(:,l);
    
    % Estimate the model without social interaction to derive one initial
    % guess
    UseX         = [ones(1,n*G);UseXc';UseXp'];
    if abs(det(UseX*(UseX')))>1.0000e-12
        betaini  = (UseX*(UseX'))\(UseX*Usey_1_true);
        sigmaini = sqrt(((Usey_1_true-(UseX')*betaini)')*(Usey_1_true-(UseX')*betaini)/(n*G-3));
        thetaini_1 = [betaini;0;sigmaini];  
    else
        thetaini_1 = [0.1 0.9 1.2 0 1.3]';
    end
    
    recordthetaini_1(:,l) = thetaini_1;
    
    % estimation under true information structure
    n0            = n;
    G0            = G;
    n1            = floor(n*n1pro);
    n10           = n1;
    FXPiter0      = iter1;
    Tol0          = TolX1;
    options       = optimset('MaxIter', iter3, 'TolFun', TolF, 'TolX', TolX3);
    theta2select  = zeros(5,4);
    logL2select   = zeros(4,1);
    exit2select   = 8*ones(4,1);
    [theta2select(:,1),logL2select(1),exit2select(1)]=fminsearch(@(theta) AvgLogLIterInitPublicDiscrete1...
        (theta,UseDDW,UseXc,UseXp,Usey_1_true,n0,G0,FXPiter0,Tol0),thetaini_1,options);
    [theta2select(:,2),logL2select(2),exit2select(2)]=fminsearch(@(theta) AvgLogLIterInitPublicDiscrete1...
        (theta,UseDDW,UseXc,UseXp,Usey_1_true,n0,G0,FXPiter0,Tol0),thetaini_2,options);
    [theta2select(:,3),logL2select(3),exit2select(3)]=fminsearch(@(theta) AvgLogLIterInitPublicDiscrete1...
        (theta,UseDDW,UseXc,UseXp,Usey_1_true,n0,G0,FXPiter0,Tol0),thetaini_3,options);
    [theta2select(:,4),logL2select(4),exit2select(4)]=fminsearch(@(theta) AvgLogLIterInitPublicDiscrete1...
        (theta,UseDDW,UseXc,UseXp,Usey_1_true,n0,G0,FXPiter0,Tol0),thetaini_4,options);
    [sortlogL,ranklogL]=sort(logL2select,'ascend');
    thetaEst(:,l)=theta2select(:,ranklogL(1)); 
    logLEst(l)=-sortlogL(1);   
    exitcond(l)=exit2select(ranklogL(1));
    select(l)=ranklogL(1);
    
    
    countloop=countloop+1;                                                 % count the number of loops calculated
    
    save MonteContinuousPublicDiscreteCompareEst.mat L G n n1 Xc Xp Wrecord yrecord beta_true lambda_true sigmasq_true...
    recordthetaini_1 thetaini_2 thetaini_3 thetaini_4 thetaEst logLEst exitcond select countloop   
end

fintime = cputime;
elapsed = toc;
time2   = clock;
fprintf('TIC TOC: %g\n', elapsed);
fprintf('CPUTIME: %g\n', fintime - initime);
fprintf('CLOCK:   %g\n', etime(time2, time1));    
Ctime=fintime-initime; 
save MonteContinuousPublicDiscreteCompareEst1.mat L G n n1pro Xc Xp yrecord beta_true lambda_true sigmasq_true...
    recordthetaini_1 thetaini_2 thetaini_3 thetaini_4 thetaEst logLEst exitcond select countloop...
    elapsed Ctime friendnum iter1 iter2 TolX1 TolX2 iter3 TolF TolX3                                                                                                                                                                                                                                                                                                                 




  

















