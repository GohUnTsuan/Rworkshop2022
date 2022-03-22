% This file is made to generate the data for the Monte Carlo experiments
% for simultaneous SAR Tobit model

% The codes about the simultaneous spatial Tobit model were used for the paper:  

% Yang, Chao, Lung-fei Lee, and Xi Qu. "Tobit models with social interactions: 
% Complete vs incomplete information."?Regional Science and Urban Economics?73 (2018): 30-50.


% Simulate the sample
clear             
tic;
initime=cputime;          
time1=clock;
rng(20160415);                         
G             = 10;                                                       % Number of independent groups
n             = 20;                                                        % Homogeneous size of each group  
beta_true     = [0 1 1];                                                   % true values of the parameters
lambda_true   = 0.6;                                                       % intensity of social interactions
sigmasq_true  = 1;                                                         % variance of the idiosyncratic risks
mu_true       = 1;                                                         % mean of privately known characteristics
etasq_true    = 4;                                                         % variance of the privately known characteristics
rho_true      = 0.4;                                                       % coefficient of correlation
Sigma_true    = etasq_true*(rho_true*ones(n)+(1-rho_true)*eye(n));         % variance-covariance matrix
R             = chol(Sigma_true);                                          
friendnum     = 15;                                                        % set the number of friends an agent can make in a grou   
iter2         = 200;
TolX2         = 1.0000e-6;                                        

% Case 1 continuous model with public information
% theta=[beta';log((1+lambda)/(1-lambda));sigma]';
% sigma=exp(theta(5));
% sigma^2=exp(2*theta(5));
% Since |lambda|<1, we make the transformation that   
% lambda=(exp(theta(4))-1)/(exp(theta(4))+1).
% So theta(4)=log((1+lambda)/(1-lambda)). 

L = 100;   
a = 1;     
for i=1:a 
    Xc            = randn(G*L,n);                                              % true values for the commonly known personal characteristics
    Xc            = reshape(Xc',n*G,L);
    Xp            = mu_true*ones(G*L,n)+randn(G*L,n)*R;                    % true values for the privately known personal characteristics
    Xp            = reshape(Xp',n*G,L);
    epsilon_base  = randn(n*G,L); 
    subW          = zeros(n*n*G,L);
    subW          = sparse(subW);
    y_3_true      = zeros(n*G,L);  
    for l=1:L  
        Wlatent         = rand(n*G,n).*(ones(n*G,n)-repmat(eye(n),G,1));       % Latent variable for the interaction structure; the self-relation latent value is set to be zero
        [Wlatentr,rank] = sort(Wlatent,2,'descend');                           % rank in descendent order
        [rankr,W]       = sort(rank,2,'ascend');                               % generate the network connections
        W               = (W<=friendnum);                                      % Everyone nominates two friends in the group  
        subW            = subW+sparse(1:(n*n*G),l*ones(1,n*n*G),reshape(W',n*n*G,1)',n*n*G,L);    % record network structure
        UseDDW          = sparse(reshape(repmat(1:(n*G),n,1),1,n*n*G),reshape(repmat(reshape(1:(n*G),n,G),n,1),1,n*n*G),...
            reshape(W',n*n*G,1),n*G,n*G);
        UseDDW          = UseDDW/friendnum; 
        UseXc           = Xc(:,l);
        UseXp           = Xp(:,l);  
        Useepsilon_base = epsilon_base(:,l);
        
        ystar_3_ini     = beta_true(1)+UseXc*beta_true(2)+beta_true(3)*UseXp;
 
        % generate the true value of the observations in that case
        
        simy_3_true     = FXP_complete_3([beta_true';log((1+lambda_true)/(1-lambda_true));(1/2)*log(sigmasq_true)],...
            UseXc,UseXp,UseDDW,Useepsilon_base,ystar_3_ini,iter2,TolX2);
        
        y_3_true(:,l)   = simy_3_true;    
    end
    filename            = sprintf('%s_%d.mat','MCContinuousPrivateOneTrue',i);     
    save(filename,'a', 'L', 'n', 'G', 'Xc', 'Xp', 'subW', 'y_3_true',...
        'beta_true', 'lambda_true', 'sigmasq_true', 'mu_true', 'etasq_true', 'rho_true', 'friendnum')   
end

fintime = cputime;
elapsed = toc;
time2   = clock;
fprintf('TIC TOC: %g\n', elapsed);
fprintf('CPUTIME: %g\n', fintime - initime);
fprintf('CLOCK:   %g\n', etime(time2, time1));      
Ctime=fintime-initime;      

                                                                                                                                                                                                                                        




  

















