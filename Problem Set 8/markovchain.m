function [con,s] = markovchain(param,c,fspace, T)
%markovchain Generate Markov chain
%   [chain] = markovchain(param,start)
s = ones(T,2);  % Initialize Markov Chain
con = ones(T,1);  % Initialize Markov Chain

s(1,:)= [0, param.ygrid(3)];    % Set starting value

cum_prob = cumsum(param.yPP,2);  % Compute cumulative distribution

% Generate Markov Chain using random numbers uniformly distributed
for t = 2:T
    con(t-1) = funeval(c,fspace,s(t-1,:));
    s(t,1) = (1+param.r)*s(t-1,1)+exp(s(t-1,2))-con(t-1);
    s(t,2) = param.ygrid(find(cum_prob(param.ygrid==s(t-1,2),:)>rand(),1));
end

end