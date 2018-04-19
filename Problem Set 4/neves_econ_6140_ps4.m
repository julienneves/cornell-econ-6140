%% Value Function Iteration + Shooting Algorithm
% Course: ECON 6140
% Version: 1.0
% Author: Julien Neves
clear, clc;

%% Question 4
n = 1000;	% Size of grid

alpha = 0.33;   % labor share
delta = 0.05;   % depreciation of capital
sigma = 0.5;    % CRRA
beta = 0.98;     % discount factor
tau_x = 0.01;   % technology
tau_c = 0.01;   % technology

crit = 1;	% Initialize convergence criterion
tol = .001;	% Convergence tolerance

T = 60;

% Grid
kstar = ((alpha*beta)/((1-beta*(1-delta))*(1+tau_x)))^(1/(1-alpha)); % Steady state
cstar = (kstar^alpha-delta*kstar*(1+tau_x))/(1+tau_c);

ub_k = 0.9*kstar;	% Upper bound
lb_k = 1.1*kstar;	% Lower bound
kgrid = linspace(lb_k,ub_k,n)';	% Create grid

k0 = 0.9*kstar;

% Empty
val_temp = zeros(n,1);  % Initialize temporary value function vector
val_fun = zeros(n,1);	% Initialize value function vector
pol_fun = zeros(n,1);	% Initialize policy function vector

ite = 0;    % Initialize iteration counter

% Value function iteration
while crit>tol ;
    % Iterate on k
    for i=1:n   
        c = (kgrid(i)^alpha + (1+tau_x)*((1-delta)*kgrid(i) - kgrid))/(1+tau_c); % Compute consumption for kt
        utility_c = c.^(1-sigma)/(1-sigma); % Compute utility for every ct
        utility_c(c<=0) = -Inf;	% Set utility to -Inf for c<=0
        [val_fun(i),pol_fun(i)] = max(utility_c + beta*val_temp);   % Solve bellman equation
    end
    crit = max(abs(val_fun-val_temp));  % Compute convergence criterion
    val_temp = val_fun;	% Update value function
    ite = ite + 1 % update iteration
end

% Value function
figure(1)
plot(kgrid,val_fun) 
xlabel('k')
ylabel('v(k)')
print -depsc fig1.eps

% Policy function
figure(2)
plot(kgrid,kgrid(pol_fun))
axis equal
hold on
plot(kgrid,kgrid, '--b')
xlabel('k')
ylabel('g(k)')
legend('Policy Function','45 degree', 'Location', 'best')
print -depsc fig2.eps


%% Question 5
kpath=zeros(1,T);
[~,kpath(1)]=min(abs(kgrid-k0));

for i = 1:T
    kpath(i+1) = pol_fun(kpath(i));
end
kpath = kgrid(kpath);

% Capital path
figure(3)
plot(0:length(kpath)-1,kpath, '--o');
xlabel('t') % x-axis label
ylabel('k') % y-axis label
print -depsc fig3.eps

%% Question 6
kstar_old = kstar;
cstar_old = cstar;

N = 100000; % grid size
T = 100;    % time periods
tol = .01;	% Convergence tolerance

tau_c = repmat(.03,1,T);   % technology
tau_x = repmat(.05,1,T);   % technology

kstar = ((alpha*beta)/((1-beta*(1-delta))*(1+tau_x(T))))^(1/(1-alpha)); % Steady state
cstar = (kstar^alpha-delta*kstar*(1+tau_x(T)))/(1+tau_c(T));   % c steady state 

k0 = kstar_old; % starting k value

k = zeros(1,T+1);   % initial k path vector
c = zeros(1,T+1);   % initial c path vector

lb_k = 0.9*kstar;   % lower bound of k axis
ub_k = 1.1*kstar;   % upper bound of k axis

lb_c = 0.9*cstar;   % lower bound of c axis
ub_c = 1.1*cstar;   % lower bound of c axis

kgrid = linspace(lb_k,ub_k,n)';   % k axis
cgrid = linspace(lb_c,ub_c,n)';   % c axis

crit = 1;   % initialize tolerance criteria
ite = 1;    % initialize iteration

while (crit>tol && ite<=N)
    k(1) = k0;  % set starting k0
    c(1) = cgrid(ite); % pick c0
    for t = 1:T-1
        k(t+1) = k(t)^alpha/(1+tau_x(t))+(1-delta)*k(t)-c(t)*(1+tau_c(t))/(1+tau_x(t)); % compute k(t+1)
        A = (1+tau_c(t))*(1+tau_x(t+1))/((1+tau_x(t))*(1+tau_c(t+1))); % tax scaling factor
        c(t+1) = c(t)*(A*(alpha*beta*k(t)^(alpha-1)/(1+tau_x(t+1))+beta*(1-delta)))^(1/sigma); %compute c(t+1)
        crit = max(abs(kstar-k(t+1)),abs(cstar-c(t+1)));    % deviation from steady state
        if crit<=tol
            % if close to steady state stop algorithm
            k = k(1:t+1); % cut path after convergences
            c = c(1:t+1); % cut path after convergences
            break
        else
            continue
        end
    end
    ite = ite + 1 % update iteration
end

k = [kstar_old, k]; % add k0
c = [cstar_old, c]; % add c0

% plot time path of k and c
figure(4)
subplot(2,1,1)
plot(0:length(k)-1,k, '--o');
xlabel('t') % x-axis label
ylabel('k') % y-axis label
subplot(2,1,2)
plot(0:length(c)-1,c, '--o');
xlabel('t') % x-axis label
ylabel('c') % y-axis label
print -depsc fig4.eps

%% Question 7
T= 500
tau_c = [repmat(.03,1,50),repmat(.01,1,T-10)];   % technology
tau_x = [repmat(.05,1,50),repmat(.01,1,T-10)];   % technology

kstar = ((alpha*beta)/((1-beta*(1-delta))*(1+tau_x(T))))^(1/(1-alpha)); % Steady state
cstar = (kstar^alpha-delta*kstar*(1+tau_x(T)))/(1+tau_c(T));   % c steady state 

k0 = kstar_old; % starting k value

k = zeros(1,T+1);   % initial k path vector
c = zeros(1,T+1);   % initial c path vector

lb_k = 0.9*kstar;   % lower bound of k axis
ub_k = 1.1*kstar;   % upper bound of k axis

lb_c = 0.9*cstar;   % lower bound of c axis
ub_c = 1.1*cstar;   % lower bound of c axis

kgrid = linspace(lb_k,ub_k,n)';   % k axis
cgrid = linspace(lb_c,ub_c,n)';   % c axis

crit = 1;   % initialize tolerance criteria
ite = 1;    % initialize iteration

while (crit>tol && ite<=N)
    k(1) = k0;  % set starting k0
    c(1) = cgrid(ite); % pick c0
    for t = 1:T-1
        k(t+1) = k(t)^alpha/(1+tau_x(t))+(1-delta)*k(t)-c(t)*(1+tau_c(t))/(1+tau_x(t)); % compute k(t+1)
        %A = (1+tau_c(t))*(1+tau_x(t+1))/((1+tau_x(t))*(1+tau_c(t+1))); % tax scaling factor
       A = 1;
        c(t+1) = c(t)*(A*(alpha*beta*k(t)^(alpha-1)/(1+tau_x(t+1))+beta*(1-delta)))^(1/sigma); %compute c(t+1)
        crit = max(abs(kstar-k(t+1)),abs(cstar-c(t+1)));    % deviation from steady state
        if crit<=tol
            % if close to steady state stop algorithm
            k = k(1:t+1); % cut path after convergences
            c = c(1:t+1); % cut path after convergences
            break
        else
            continue
        end
    end
    ite = ite + 1 % update iteration
end

k = [kstar_old, k]; % add k0
c = [cstar_old, c]; % add c0

% plot time path of k and c
figure(5)
subplot(2,1,1)
plot(0:length(k)-1,k, '--o');
xlabel('t') % x-axis label
ylabel('k') % y-axis label
subplot(2,1,2)
plot(0:length(c)-1,c, '--o');
xlabel('t') % x-axis label
ylabel('c') % y-axis label
print -depsc fig5.eps
