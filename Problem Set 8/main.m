clear; clc;
%% Problem 1
% declare parameters
param.ro=0.90;
param.se=sqrt(0.06);

% part a) - b)
param.k = 5;
[param, ro_ta5, se_ta5] = transition(param, "tauchen");
param.k = 10;
[param, ro_ta10, se_ta10] = transition(param,  "tauchen");

tab_1b = table([param.ro; ro_ta5; ro_ta10],[param.se; se_ta5; se_ta10],'VariableNames',{'rho','se'},'RowNames',{'Model';'Five';'Ten'})

% part c)
param.k = 5;
[param, ro_rw, se_rw] = transition(param,  "rouwenhorst");

tab_1c = table([param.ro; ro_ta5; ro_rw],[param.se; se_ta5; se_rw],'VariableNames',{'rho','se'},'RowNames',{'Model'; 'tauchen'; 'rouwenhorst'})

% part d)
old_var = param.se^2/(1-param.ro^2);
param.ro = 0.98;
param.se = old_var*(1-param.ro^2);
[param, ro_rw, se_rw] = transition(param,  "rouwenhorst");

tab_1d = table([param.ro;ro_rw],[param.se;se_rw],'VariableNames',{'rho','se'},'RowNames',{'Model'; 'rouwenhorst'})

%% Problem 2
% declare parameters
param.beta=0.95;
param.r=0.02;
param.ro=0.9;
param.se=sqrt(0.06);
param.gamma=2;
param.amin = 0;
param.n = 25;
param.k = 5;

% part a)
[param] = transition(param,  "rouwenhorst");
[param,c,fspace,s,smin,smax] = policy_ip(param);

close all
sfine=gridmake(nodeunif(param.n*2,smin(1),smax(1)),param.ygrid);
xfine=funeval(c,fspace,sfine);

figure(1)
subplot(2,1,1)
sfine=gridmake(nodeunif(param.n*4,smin(1),smax(1)),0); %ygrid(floor(k/2)+2));
xfine=funeval(c,fspace,sfine);
plot(sfine(:,1),xfine)
xlabel({'$a$'},'Interpreter','latex')
ylabel({'$c(a,\bar{y})$'},'Interpreter','latex')
title({'Consumption policy function, $y=\bar{y}$'},'Interpreter','latex')
set(gca,'FontSize',8);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(param.k*4,smin(2),smax(2)));
xfine=funeval(c,fspace,sfine);
plot(exp(sfine(:,2)),xfine)
xlabel('$y$','Interpreter','latex')
ylabel('$c(0,y)$','Interpreter','latex')
title({'Consumption policy function, $a=0$'},'Interpreter','latex')
set(gca,'FontSize',8);
print -depsc fig1.eps

figure(2)
subplot(2,1,1)
sfine=gridmake(nodeunif(param.k*4,smin(1),smax(1)),0);
xfine=funeval(c,fspace,sfine);
plot(sfine(:,1),(1+param.r)*sfine(:,1)+exp(sfine(:,2))-xfine)
xlabel('$a$','Interpreter','latex')
ylabel('$a^{\prime}(a,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $y=\bar{y}$'},'Interpreter','latex')
set(gca,'FontSize',8);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(param.k*4,smin(2),smax(2)));
xfine=funeval(c,fspace,sfine);
plot(exp(sfine(:,2)),(1+param.r)*sfine(:,1)+exp(sfine(:,2))-xfine)
xlabel('$y$','Interpreter','latex')
ylabel('$a^{\prime}(0,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $a=0$'},'Interpreter','latex')
set(gca,'FontSize',8);
print -depsc fig2.eps

% part b)
gamma = [1,2,5];
for i = 1:3
    param.gamma = gamma(i);
    
    [param] = transition(param,  "rouwenhorst");
    [param,c,fspace] = policy_ip(param);
    
    con = markovchain(param,c,fspace, 10000); % Generate Markov chain
    se_c(i) = std(con);
end

tab_2b = table(se_c','VariableNames',{'std_c'},'RowNames',{'Gamma = 1';'Gamma = 2';'Gamma = 5'})


% part c)
param.gamma = 2;
se = [0.01,0.06,0.12];

figure(3)
hold on
for i = 1:3
    param.se = sqrt(se(i));
    
    [param] = transition(param,  "rouwenhorst");
    [param,c,fspace] = policy_ip(param);
    
    sfine=gridmake(0,nodeunif(param.k*4,smin(2),smax(2)));
    xfine=funeval(c,fspace,sfine);
    plot(exp(sfine(:,2)),1-xfine./exp(sfine(:,2)))
end
xlabel('$y$','Interpreter','latex')
ylabel('$a^{\prime}(0,y)/y$','Interpreter','latex')
title({'Savings rate, $a=0$'},'Interpreter','latex')
legend('\sigma_e^2 = 0.01','\sigma_e^2 = 0.06','\sigma_e^2 = 0.12')
set(gca,'FontSize',8);
hold off
print -depsc fig3.eps


% part d) - e)
% no-borrowing
param.gamma = 2;
param.se = sqrt(0.06);
param.amin = 0;

[param] = transition(param,  "rouwenhorst");
[param,c,fspace] = policy_ip(param);

[con, s] = markovchain(param,c,fspace, 10000);
avg_c(1) = mean(con);

c_t = log(con(2:end))-log(con(1:end-1));
e_t = s(2:end,2)-param.ro*s(1:end-1,2);
sig = cov(c_t,e_t);
phi(1) = 1-sig(1,2)/sig(2,2);

% natural debt limit
param.amin = -min(exp(param.ygrid)+.01)/param.r;

[param] = transition(param,  "rouwenhorst");
[param,c,fspace] = policy_ip(param);

[con, s] = markovchain(param,c,fspace, 10000);
avg_c(2) = mean(con);

c_t = log(con(2:end))-log(con(1:end-1));
e_t = s(2:end,2)-param.ro*s(1:end-1,2);
sig = cov(c_t,e_t);
phi(2) = 1-sig(1,2)/sig(2,2);

tab_2d = table(avg_c',phi','VariableNames',{'mean_c','phi'},'RowNames',{'No borrowing';'Natural'})


%% Problem 3
% declare parameters
param.beta=0.95;
param.r=0.02;
param.ro=0;
param.se=sqrt(0.06);
param.gamma=2;
param.amin = 0;
param.n = 25;
param.k = 7;

% part b)
wbar = -param.se^2/2;
[x,w] = qnwnorm(7,wbar,param.se^2);
fprintf('E(y)=%f \n',exp(x)'*w)

% part c)
% [param] = transition(param,  "rouwenhorst");
param.ws = w';
param.ygrid = x;
[param,c,fspace,s,smin,smax] = policy_ca(param);

figure(5)
sfine=gridmake(nodeunif(param.k*4,smin(1),smax(1)));
xfine=funeval(c,fspace,sfine);
plot(sfine,xfine)
xlabel('$x$','Interpreter','latex')
ylabel('$c(x)$','Interpreter','latex')
title({'Consumption policy function'},'Interpreter','latex')
set(gca,'FontSize',8);
print -depsc fig4.eps