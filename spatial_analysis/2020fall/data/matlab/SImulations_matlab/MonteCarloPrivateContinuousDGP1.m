% This file is made for the Monte Carlo experiment of the model of social
% interactions with private information. Consider the case that some
% characteristics are known only to an agent herself and that
% characteristics are discretely distributed.  

% The quadrature method is used to derive (numerically) the equilibrium
% conditional expectation functions.

% Codes are updated by Chao Yang, June 2019.

% REFERENCE:
% Yang and Lee, 2017, "Social interactions under incomplete information with
% heterogeneous expectations."?Journal of Econometrics?198.1 (2017): 65-83

% Simulate the sample
clear             
tic;
initime=cputime;    
time1=clock;
rng(20130622);                         
G             = 5;                                                       % Number of independent groups
n             = 20;                                                        % Homogeneous size of each group  
K             = 6; 
[quadx,quadw] = lgwt(K,-1,1);                                              % derive quadrature points and weights
ExpM0         = zeros(K*n*G,1);
beta_true     = [0 1 1];                                                   % true values of the parameters
lambda_true   = 0.3;                                                       % intensity of social interactions
sigmasq_true  = 1;                                                         % variance of the idiosyncratic risks
mu_true       = 1;                                                         % mean of privately known characteristics
etasq_true    = 4;                                                         % variance of the privately known characteristics
rho_true      = 0.4;                                                       % coefficient of correlation
Sigma_true    = etasq_true*(rho_true*ones(n)+(1-rho_true)*eye(n));         % variance-covariance matrix
R             = chol(Sigma_true);                                          
friendnum     = 3;                                                         % set the number of friends an agent can make in a grou   
iter2         = 200;
TolX2         = 1.0000e-6;                                          

% Case 1 continuous model with public information
% theta=[beta';log((1+lambda)/(1-lambda));sigma]';
% sigma=abs(theta(5));
% sigma^2=(theta(5))^2;
% Since |lambda|<1, we make the transformation that   
% lambda=(exp(theta(4))-1)/(exp(theta(4))+1).
% So theta(4)=log((1+lambda)/(1-lambda)). 

L = 20;   
a = 1;   
for i=1:a 
    Xc            = randn(G*L,n);                                              % true values for the commonly known personal characteristics
    Xc            = reshape(Xc',n*G,L);
    Xp            = mu_true*ones(G*L,n)+randn(G*L,n)*R;                    % true values for the privately known personal characteristics
    Xp            = reshape(Xp',n*G,L);
    epsilon_true  = (sqrt(sigmasq_true))*randn(n*G,L); 
    subW          = zeros(n*n*G,L);
    subW          = sparse(subW);
    ExpMrecord    = zeros(n*G,L);
    y_1_true      = zeros(n*G,L);  
    for l=1:L  
        Wlatent         = rand(n*G,n).*(ones(n*G,n)-repmat(eye(n),G,1));       % Latent variable for the interaction structure; the self-relation latent value is set to be zero
        [Wlatentr,rank] = sort(Wlatent,2,'descend');                           % rank in descendent order
        [rankr,W]       = sort(rank,2,'ascend');                               % generate the network connections
        W               = (W<=friendnum);                                      % Everyone nominates two friends in the group  
        subW            = subW+sparse(1:(n*n*G),l*ones(1,n*n*G),reshape(W',n*n*G,1)',n*n*G,L);    % record network structure
        UseDDW          = sparse(reshape(repmat(1:(n*G),n,1),1,n*n*G),reshape(repmat(reshape(1:(n*G),n,G),n,1),1,n*n*G),...
            reshape(W',n*n*G,1),n*G,n*G);
        UseDDW          = UseDDW/friendnum; 
        UseDW           = sparse(reshape(repmat(1:(n*K*G),n,1),1,n*n*K*G),reshape(repmat(reshape(1:(n*K*G),n,K*G),n,1),1,n*n*K*G),...
            reshape(repmat(reshape(W',n*n,G),K,1),1,n*n*K*G),n*K*G,n*K*G);
        UseDW           = UseDW/friendnum;
        UseRDW          = sparse(reshape(repmat(1:(n*K*G),n,1),1,n*n*K*G),reshape(repmat(reshape(1:(n*K*G),n*K,G),n,1),1,n*n*K*G),...
            reshape(repmat(reshape(W',n,n*G),K,1),1,n*K*n*G),n*K*G,n*K*G);
        UseRDW          = UseRDW/friendnum;
        UseXc           = Xc(:,l);
        UseXp           = Xp(:,l);  
        Useepsilon_true = epsilon_true(:,l);
        
        % calculate the equilibrium expectation in the true parameter values
        ExpM_1_true     = Expect_private_data_1([beta_true';log((1+lambda_true)/(1-lambda_true))],...
            mu_true,etasq_true,rho_true,UseDW,UseRDW,UseXc,UseXp,n,G,K,quadx,quadw,iter2,TolX2,ExpM0);
        
        % generate the true value of the observations in that case
        simy_1_true     = beta_true(1)+UseXc*beta_true(2)+UseXp*beta_true(3)+...
            lambda_true*ExpM_1_true-Useepsilon_true;     
        ExpMrecord(:,l) = ExpM_1_true; 
        y_1_true(:,l)   = simy_1_true;  
    end
    filename            = sprintf('%s_%d.mat','MonteLinearPrivateContinuousDGP1',i);     
    save(filename,'a', 'L', 'n', 'G', 'K','Xc', 'Xp', 'subW', 'y_1_true', 'ExpMrecord',...
        'beta_true', 'lambda_true', 'sigmasq_true', 'mu_true', 'etasq_true', 'rho_true', 'friendnum')   
end

fintime = cputime;
elapsed = toc;
time2   = clock;
fprintf('TIC TOC: %g\n', elapsed);
fprintf('CPUTIME: %g\n', fintime - initime);
fprintf('CLOCK:   %g\n', etime(time2, time1));          
Ctime=fintime-initime;                                     
                                                                                                                                                                                                                                        




  

















