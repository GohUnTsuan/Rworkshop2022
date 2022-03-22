% This m-file is written to define the total likelihood function for the
% MLE estimation of the continous choice model wirh private information and
% discrete characteristics

function [logL]=AvgLogLIterInitPublicDiscrete1(theta,UseDDW,UseXc,UseXp,Usey_1_true,n0,G0,FXPiter0,Tol0)

% calculate expectation by contraction mapping iteration
ExpMest=FXP_public_discrete_Matrix_1(theta, UseDDW, UseXc, UseXp, FXPiter0, Tol0);

% calculate (negative) sample average log likelihood
meanPrivate1=theta(1)+UseXc*theta(2)+UseXp*theta(3)+...
        ((exp(theta(4))-1)/(exp(theta(4))+1))*UseDDW*ExpMest;
logL=-(-n0*G0*log(theta(5)^2)/2-n0*G0*(1/2)*log(2*pi)-(1/(2*(theta(5)^2)))*...
    ((meanPrivate1-Usey_1_true)')*(meanPrivate1-Usey_1_true));                     








    