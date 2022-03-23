% This m-file is written to define the total likelihood function for the
% MLE estimation of the continous choice model wirh private information and
% continuous characteristics

function [logL]=continuousPrivatecondlogL1(theta,mu,etasq,rho,UseDW,UseRDW,...  
    UseXc,UseXp,Usey_1_true,n0,G0,K0,quadx,quadw,FXPiter0,Tol0) 

% calculate expectation by contraction mapping iteration
ExpMest    = FXP_private_Matrix_1(theta,mu,etasq,rho,UseDW,UseXc,n0,G0,K0,quadx,quadw,FXPiter0,Tol0);

% calculate (negative) sample average log likelihood
outcome    = theta(1)+reshape(repmat(reshape(UseXc,n0,G0),K0,1),n0*K0*G0,1)*theta(2)...
        +theta(3)*repmat(reshape(repmat(log((quadx+1)./(1-quadx))',n0,1),n0*K0,1),G0,1)+...
        ((exp(theta(4))-1)/(exp(theta(4))+1))*UseDW*ExpMest;

adjcoeff   = quadw.*exp(-((log((quadx+1)./(1-quadx))-(1-rho)*mu).^2)/(2*(1-rho^2)*etasq))...
        ./((1-quadx).*(1+quadx)); 
    
selfweight = exp(-rho*rho*(UseXp.^2)/(2*(1-rho^2)*etasq));
    
crossterm  = exp(rho*UseXp*((log((quadx+1)./(1-quadx))-(1-rho)*mu)')/((1-rho^2)*etasq));
    
adjweight  = sparse(1:(n0*G0),1:(n0*G0),selfweight,n0*G0,n0*G0)*crossterm*sparse(1:K0,1:K0,adjcoeff,K0,K0); 

ExpMdata   =  sqrt(2/(pi*(1-rho^2)*etasq))*sparse(reshape(repmat(1:(n0*G0),K0,1),1,n0*K0*G0),...
    1:(n0*K0*G0),ones(1,n0*K0*G0),n0*G0,n0*K0*G0)*sparse(1:(n0*K0*G0),1:(n0*K0*G0),...
    reshape(adjweight',1,n0*K0*G0),n0*K0*G0,n0*K0*G0)*UseRDW*outcome;

meanPrivate1 = theta(1)+UseXc*theta(2)+UseXp*theta(3)+...
        ((exp(theta(4))-1)/(exp(theta(4))+1))*ExpMdata;

logL=-(-n0*G0*log(theta(5)^2)/2-n0*G0*(1/2)*log(2*pi)-(1/(2*(theta(5)^2)))*...
    ((meanPrivate1-Usey_1_true)')*(meanPrivate1-Usey_1_true));                                   








    