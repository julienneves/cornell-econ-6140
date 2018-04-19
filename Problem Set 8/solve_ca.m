function [x] = solve_ca(param,c, fspace, s)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

ns=length(s);

a=.01*ones(ns,1);
b=(1+param.r)*s(:,1)+exp(s(:,2));
tol=1e-8; %tolerance level

fa=euler_ip(a,c,fspace,s,param);
fb=euler_ip(b,c,fspace,s,param);

x=zeros(ns,1);

% Start bisection
dx = 0.5*(b - a);

x = a + dx;                       %  start at midpoint
sb=sign(fb);
dx = sb.*dx;                      %  we increase or decrease x depending if f(b) is positive or negative                     

i=0;
  while any(abs(dx)>tol)
   i=i+1;
    dx = 0.5*dx;
    x = x - sign(euler_ip(x,c,fspace,s,param)).*dx;   
  end


x(fb>=0)=b(fb>=0);
end

