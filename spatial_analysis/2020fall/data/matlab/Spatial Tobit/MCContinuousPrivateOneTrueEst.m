
% This file is made for the Monte Carlo experiment 

% First, load the sample       
clear                            
tic;
initime          = cputime;             
time1            = clock;
rng(20160424);
load('MCContinuousPrivateOneTrue_1.mat')                                                                      
iter3           = 2500;                   
TolF            = 1.0000e-6;
TolX3           = 1.0000e-6;         


% reparametrization
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

% estimate the data under models of different information structures

thetaXpest       = zeros(3,L);

recordthetaini_1 = zeros(5,L);                                             % store the initial guess by regression without network effects
thetaini_2       = 0.5*rand(5,1);                                          % initial guess randomly chosen
thetaini_3       = [0 0 0 0 1]';                                           % initial guess using 0 and 1  

thetaEstcom      = repmat([0 0 0 0.5 0.8]',1,L); 
logLEstcom       = ones(L,1);
exitcondcom      = 8*ones(L,1);
selectcom        = zeros(L,1); 
countloop_1      = 0; 
for l=1:L 
    
    W            = reshape(subW(:,l),n,n*G)';
    UseDDW       = sparse(reshape(repmat(1:(n*G),n,1),1,n*n*G),reshape(repmat(reshape(1:(n*G),n,G),n,1),1,n*n*G),...
         reshape(W',n*n*G,1),n*G,n*G);
    UseDDW       = UseDDW/friendnum; 
    UseXc        = Xc(:,l);
    UseXp        = Xp(:,l);                                                         % select sample used in the l-th repetition
    Usey_3_true  = y_3_true(:,l);
    options       = optimset('MaxIter', iter3, 'TolFun', TolF, 'TolX', TolX3);  
    
    
    % Estimate the model without social interaction to derive one initial
    % guess
    UseX             = [ones(1,n*G);UseXc';UseXp']';
    fsthetaEst       = fminunc(@(fstheta) fslogLtobit(fstheta,UseX,Usey_3_true),[0 0 0 1]',options);
    thetaini_1       = [fsthetaEst(1:3);0;fsthetaEst(4)];  
     
    recordthetaini_1(:,l) = thetaini_1;   
    
    % estimation under complete information
    
    theta2selectcom   = zeros(5,3);
    logL2selectcom    = zeros(3,1);
    exit2selectcom    = 8*ones(3,1);
    [theta2selectcom(:,1),logL2selectcom(1),exit2selectcom(1)] = fminsearch(@(theta) continuousCompletecondlogL3...
        (theta,UseXc,UseXp,Usey_3_true,UseDDW),thetaini_1,options);
    [theta2selectcom(:,2),logL2selectcom(2),exit2selectcom(2)] = fminsearch(@(theta) continuousCompletecondlogL3...
        (theta,UseXc,UseXp,Usey_3_true,UseDDW),thetaini_2,options);
    [theta2selectcom(:,3),logL2selectcom(3),exit2selectcom(3)] = fminsearch(@(theta) continuousCompletecondlogL3...
        (theta,UseXc,UseXp,Usey_3_true,UseDDW),thetaini_3,options);
    [sortlogLcom,ranklogLcom]=sort(logL2selectcom,'ascend');
    thetaEstcom(:,l)   = theta2selectcom(:,ranklogLcom(1)); 
    logLEstcom(l)      = -sortlogLcom(1);   
    exitcondcom(l)     = exit2selectcom(ranklogLcom(1));
    selectcom(l)       = ranklogLcom(1);
    
    
    
    
    countloop_1       = countloop_1+1;                                                 % count the number of loops calculated
    
    

    
end



fintime = cputime;
elapsed = toc;
time2   = clock;
fprintf('TIC TOC: %g\n', elapsed);
fprintf('CPUTIME: %g\n', fintime - initime);
fprintf('CLOCK:   %g\n', etime(time2, time1));   
Ctime=fintime-initime; 
save MCContinuousPrivateOneTrueEst_1.mat L G n Xc Xp subW y_3_true...
    beta_true lambda_true sigmasq_true mu_true etasq_true rho_true...
    recordthetaini_1 thetaini_2 thetaini_3...
    thetaEstcom logLEstcom exitcondcom selectcom...
    friendnum iter3 TolF TolX3 countloop_1                    
    
    











  

















