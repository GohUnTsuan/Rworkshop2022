% This file is made for the Monte Carlo experiment of the model of social
% interactions with private information. This set of experiments consider
% the case that all the variables available in the data sets are also
% publcly known to all group members.

% Codes are updated by Chao Yang, June 2019.

% REFERENCE:
% Yang and Lee, 2017, "Social interactions under incomplete information with
% heterogeneous expectations."?Journal of Econometrics?198.1 (2017): 65-83.
% 

% Simulate the sample
clear             
tic;
initime=cputime;  
time1=clock;
rng(20130622);                         
L             = 100;                                                      % Number of simulations
G             = 10;                                                         % Number of independent groups
n             = 20;                                                      % Homogeneous size of each group  
ExpM0         = zeros(n*G,1);
beta_true     = [0 1 1];                                                   % true values of the parameters
lambda_true   = 0.3;                                                       % intensity of social interactions
sigmasq_true  = 1;                                                         % variance of the idiosyncratic risks
Xc            = randn(G*L,n);                                              % true values for the commonly known personal characteristics
Xc            = reshape(Xc',n*G,L);
n1pro         = 0.4;  
n1            = floor(n*n1pro);                                            % the number of agents who satisfy a criterion
Xplatent      = n*rand(G*L,n);                                             % generate the true value of the privately known personal characteristics 
[sortXplatent,orderXplatent]  = sort(Xplatent,2,'descend');
[sortorder,Xprank]            = sort(orderXplatent,2,'ascend');
Xp            = (Xprank<=n1); 
Xp            = reshape(Xp',n*G,L);
epsilon_true  = (sqrt(sigmasq_true))*randn(n*G,L); 
friendnum     = 3;                                                        % set the number of friends an agent can make in a grou   
iter2         = 200;
TolX2         = 1.0000e-6;                                     

% Case 1 continuous model with public information
% theta=[beta';log((1+lambda)/(1-lambda));sigma]';
% sigma=abs(theta(5));
% sigma^2=(theta(5))^2;
% Since |lambda|<1, we make the transformation that 
% lambda=(exp(theta(4))-1)/(exp(theta(4))+1).
% So theta(4)=log((1+lambda)/(1-lambda)).  

Wrecord      = zeros(n*n*G,L);
Wrecord      = sparse(Wrecord);
ExpMrecord   = zeros(n*G,L);
yrecord      = zeros(n*G,L);
countloop    = 0;   
%countloopËæÊ±´¢´æ½á¹û
for l=1:L  
    Wlatent         = rand(n*G,n).*(ones(n*G,n)-repmat(eye(n),G,1));       % Latent variable for the interaction structure; the self-relation latent value is set to be zero
    [Wlatentr,rank] = sort(Wlatent,2,'descend');                           % rank in descendent order
    [rankr,W]       = sort(rank,2,'ascend');                               % generate the network connections
    W               = (W<=friendnum);                                      % Everyone nominates two friends in the group  
    Wrecord         = Wrecord+sparse(1:(n*n*G),l*ones(1,n*n*G),reshape(W',n*n*G,1)',n*n*G,L);    % record network structure
    UseDDW          = sparse(reshape(repmat(1:(n*G),n,1),1,n*n*G),reshape(repmat(reshape(1:(n*G),n,G),n,1),1,n*n*G),...
        reshape(W',n*n*G,1),n*G,n*G);
    UseDDW          = UseDDW/friendnum; 
    UseXc           = Xc(:,l);
    UseXp           = Xp(:,l);  
    Useepsilon_true = epsilon_true(:,l);
    ExpM_1_true     = FXP_public_discrete_1([beta_true';log((1+lambda_true)/(1-lambda_true))],...
        UseDDW,UseXc,UseXp,iter2,TolX2,ExpM0);                               % calculate the equilibrium expectation in the true parameter values
    y_1_true        = beta_true(1)+UseXc*beta_true(2)+UseXp*beta_true(3)+...
        lambda_true*UseDDW*ExpM_1_true-Useepsilon_true;     % generate the true value of the observations in that case
    ExpMrecord(:,l) = ExpM_1_true; 
    yrecord(:,l)    = y_1_true;
    
    countloop       = countloop+1;                                                  % count the number of loops calculated
    
    save MonteContinuousPublicDiscreteCompare.mat L G n n1 Xc Xp Wrecord epsilon_true ExpMrecord yrecord countloop   
end

fintime = cputime;
elapsed = toc;
time2   = clock;
fprintf('TIC TOC: %g\n', elapsed);
fprintf('CPUTIME: %g\n', fintime - initime);
fprintf('CLOCK:   %g\n', etime(time2, time1)); 
Ctime=fintime-initime; 
save MonteContinuousPublicDiscreteCompare1.mat L G n n1 Xc Xp Wrecord epsilon_true ExpMrecord yrecord countloop...
    elapsed Ctime friendnum iter2 TolX2 beta_true lambda_true sigmasq_true n1pro                                                                                                                                                                                                                                         




  

















