
% This file is writtn to calculate the log likelihood function which is
% used to derive one initial guess of parameter values for Tobit model

function [fslogL] = fslogLtobit(fstheta,UseX,Usey_3_true)  
outcome           = UseX*fstheta(1:3);
p0                = 1-normcdf(outcome/exp(fstheta(4)));
p0                = p0.*(p0>1.0000e-312)+(1.0000e-312)*(p0<=1.0000e-312);
fslogL            = -(((Usey_3_true<=0)')*log(p0)+((Usey_3_true>0)')*...
    (-fstheta(4)-(1/2)*log(2*pi)-(1/(2*exp(2*fstheta(4))))*...
    ((outcome-Usey_3_true).^2)));                                                   