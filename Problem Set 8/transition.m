function [param, ro_tilde, se_tilde] = transition(param, method)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

var_w= param.se^2/(1-param.ro^2);
wbar = -var_w*(1-param.ro)/2;
% nodes for the distribution of income shocks

if method=="tauchen"
    [e,w] = tauchen(param.ro,param.se, param.k, 3);
elseif method == "rouwenhorst"
    [e,w] = rouwenhorst(param.ro,param.se,param.k);        % Rouwenhorst method
end

param.yPP=w;
param.ygrid=e';
param.ygrid = param.ygrid - wbar;

[~,D,V]=eig(w);
ws = V(:,1);
ws = ws/sum(ws);

param.ws = ws;
param.ygrid = param.ygrid - wbar;

ro_tilde = e*(w*e')/(e*e');
se_tilde = sqrt(e.^2*ws*(1-ro_tilde^2));
end

