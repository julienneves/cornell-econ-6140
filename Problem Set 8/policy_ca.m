function [param,c,fspace,s,smin,smax] = policy_ca(param)
% Bounds for state space
ymin=min(param.ygrid);
ymax=max(param.ygrid);

xmin = exp(ymin);                        % no borrowing
xmax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it

% Declare function space to approximate a'(a,y)
n=param.n;

% Lower and higher bound for the state space (a,y)
smin=xmin;
smax=xmax;

scale=1/2;
fspace=fundef({'spli', nodeunif(n(1),(.01).^scale,(xmax-xmin+.01).^scale).^(1/scale)+xmin-.01,0,3});

grid=funnode(fspace);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)

c=funfitxy(fspace,s,s/2);                %guess that keep constant assets

tic
for it=1:101
    cnew=c;
    x = solve_ca(param,c,fspace,s);
    c=funfitxy(fspace,s,x);
    
    fprintf('%4i %6.2e\n',[it,norm(c-cnew)]);
    if norm(c-cnew)<1e-5, break, end
end
toc
end

