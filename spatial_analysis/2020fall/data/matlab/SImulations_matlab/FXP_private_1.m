% This m-file is written for the fixed point iteration algorithm in the
% case of continuous choice with private information where characteristics
% are continuously distributed.

% The contraction mapping iterations are used to derive the unique fixed
% point

% Updated by Chao Yang, June 2019.

function [ExpM_1]=FXP_private_1(theta,mu,etasq,rho,VDW,VXc,n0,G0,K0,quadx,quadw,FXPiter0,Tol0,ExpM0)      
 
index=0;
i=1;
ExpM_1_0=ExpM0;  
while index<1
    outcome    = theta(1)+reshape(repmat(reshape(VXc,n0,G0),K0,1),n0*K0*G0,1)*theta(2)...
        +theta(3)*repmat(reshape(repmat(log((quadx+1)./(1-quadx))',n0,1),n0*K0,1),G0,1)+...
        ((exp(theta(4))-1)/(exp(theta(4))+1))*VDW*ExpM_1_0;  
    
    adjcoeff   = quadw.*exp(-((log((quadx+1)./(1-quadx))-(1-rho)*mu).^2)/(2*(1-rho^2)*etasq))...
        ./((1-quadx).*(1+quadx)); 
    
    selfweight = exp(-rho*rho*(log((quadx+1)./(1-quadx)).^2)/(2*(1-rho^2)*etasq));
    
    crossterm  = exp(rho*log((quadx+1)./(1-quadx))*((log((quadx+1)./(1-quadx))-(1-rho)*mu)')/((1-rho^2)*etasq));
    
    adjweight  = sparse(1:K0,1:K0,selfweight,K0,K0)*crossterm*sparse(1:K0,1:K0,adjcoeff,K0,K0); 
    
    ExpM_1_1   = sqrt(2/(pi*(1-rho^2)*etasq))*sparse(reshape(repmat(reshape(1:(n0*K0*G0),n0,K0*G0),K0,1),1,n0*K0*K0*G0),...
        reshape(repmat(reshape(1:(n0*K0*G0),n0*K0,G0),K0,1),1,n0*K0*K0*G0),...
        repmat(reshape(repmat(reshape(adjweight',1,K0*K0),n0,1),1,n0*K0*K0),1,G0),n0*K0*G0,n0*K0*G0)*outcome;  
    
    if max(abs(ExpM_1_1-ExpM_1_0))>=Tol0 && i<=FXPiter0    
        ExpM_1_0=ExpM_1_1;  
        i=i+1;
        index=0;
    else
        index=1;
    end
end
ExpM_1=ExpM_1_0;           


    