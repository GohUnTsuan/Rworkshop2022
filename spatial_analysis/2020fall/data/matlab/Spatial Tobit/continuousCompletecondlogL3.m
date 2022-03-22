% This m-file is written for the log likelihood under complete information

function [logL]=continuousCompletecondlogL3(theta,UseXc,UseXp,Usey_3,UseDDW) 

meancomp  = theta(1)+UseXc*theta(2)+UseXp*theta(3)+...
    ((exp(theta(4))-1)/(exp(theta(4))+1))*UseDDW*Usey_3; 

prob_0    = 1-normcdf(meancomp/exp(theta(5))); 
prob_0    = prob_0.*(prob_0>1.0000e-312)+(1.0000e-312)*(prob_0<=1.0000e-312); 

index     = (Usey_3 > 0);

nuncensor = sum(index);

UseDDW22  = UseDDW(index,index);
    
logL      = -((1-index')*log(prob_0)+...
    (index')*(-theta(5)-(1/2)*log(2*pi)-(1/(2*exp(2*theta(5))))*((meancomp-Usey_3).^2))+...
    log(abs(det(eye(nuncensor)-((exp(theta(4))-1)/(exp(theta(4))+1))*UseDDW22))));  

                              








    