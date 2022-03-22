

function [logLXp]=xplogL(theta2,UseXp,n0,G0)    

etasq         = theta2(2)^2;
rho           = (exp(theta2(3))-1)/(exp(theta2(3))+1);
invSigma      = (1/(1-rho))*eye(n0)-(rho/((1-rho)*(1+(n0-1)*rho)))*ones(n0);
diaginvSigma  = sparse(reshape(repmat(1:(n0*G0),n0,1),1,n0*n0*G0),...
    reshape(repmat(reshape(1:(n0*G0),n0,G0),n0,1),1,n0*n0*G0),...
    repmat(reshape(invSigma',1,n0*n0),1,G0),n0*G0,n0*G0);
logLXp        = -(-n0*G0*(1/2)*log(2*pi)-(G0/2)*log(1+(n0-1)*rho)...
    -((n0-1)*G0/2)*log(1-rho)-(n0*G0/2)*log(etasq)...
    -(1/(2*etasq))*((UseXp-theta2(1))')*diaginvSigma*(UseXp-theta2(1)));    
  
    
    