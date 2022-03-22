
% In this m file the equilibrium outcome of the simultaneous spatial tobit
% model are solved as fixed points from contraction mapping iterations 

% Codes are updated by Chao Yang.

function [y_3] = FXP_complete_3(theta,VXc,VXp,VDDW,epsilon0,ystar_3_ini,FXPiter0,Tol0)

 index             = 0;
 i                 = 1;
 ystar_3_0         = ystar_3_ini;
 
while index<1
    
    y_3_iter       = ystar_3_0.*(ystar_3_0 > 0);
    
    ystar_3_1      = theta(1)+VXc*theta(2)+theta(3)*VXp+...
        ((exp(theta(4))-1)/(exp(theta(4))+1))*VDDW*y_3_iter...
        +exp(theta(5))*epsilon0;
    
    if max(abs(ystar_3_1-ystar_3_0))>=Tol0 && i<=FXPiter0    
        ystar_3_0 = ystar_3_1;  
        i         = i +1;
        index     = 0;
    else
        index     = 1;
    end
end
ystar_3          = ystar_3_0;     
y_3              = ystar_3.*(ystar_3> 0);      

end

