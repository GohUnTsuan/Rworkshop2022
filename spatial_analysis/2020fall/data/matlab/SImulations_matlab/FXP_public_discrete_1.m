% This m-file is written for the fixed point iteration algorithm in the
% case of continuous choice with public information.

% The contraction iteration algorithm is used to derive the unique fixed
% point.

% Codes are writted and updated by Chao Yang.

function [ExpM_1]=FXP_public_discrete_1(theta, VDDW, VXc, VXp, FXPiter0, Tol0, ExpM0)
index=0;
i=1;
ExpM_1_0=ExpM0;  
while index<1
    ExpM_1_1=theta(1)+VXc*theta(2)+theta(3)*VXp+...
        ((exp(theta(4))-1)/(exp(theta(4))+1))*VDDW*ExpM_1_0;  
    if max(abs(ExpM_1_1-ExpM_1_0))>=Tol0 && i<=FXPiter0    
        ExpM_1_0=ExpM_1_1;  
        i=i+1;
        index=0;
    else
        index=1;
    end
end
ExpM_1=ExpM_1_0;                                                
    