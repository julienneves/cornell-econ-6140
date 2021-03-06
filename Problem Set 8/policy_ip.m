function [param,c,fspace,s,smin,smax] = policy_ip(param)
% Bounds for state space
ymin=min(param.ygrid);
ymax=max(param.ygrid);

amin = param.amin;                        % no borrowing
amax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it

% Declare function space to approximate a'(a,y)
n=[param.n,param.k];

% Lower and higher bound for the state space (a,y)
smin=[amin,ymin];
smax=[amax,ymax];

scale=1/2;
fspace=fundef({'spli', nodeunif(n(1),(smin(1)-amin+.01).^scale,(smax(1)-amin+.01).^scale).^(1/scale)+amin-.01,0,3},...
    {'spli', param.ygrid,0,1});

grid=funnode(fspace);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)

c=funfitxy(fspace,s,param.r/(1+param.r)*s(:,1)+exp(s(:,2)));                %guess that keep constant assets

tic
for it=1:101
    cnew=c;
    x = solve_ip(param,c,fspace,s, amin);
    c=funfitxy(fspace,s,x);
    
    fprintf('%4i %6.2e\n',[it,norm(c-cnew)]);
    if norm(c-cnew)<1e-7, break, end
end
toc
end

