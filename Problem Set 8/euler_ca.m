function fval = euler_ca(x,c,fspace,s,param)

ns=size(s,1);
xprime = (1+param.r)*(s(:,1)-x)+exp(s(:,2));
fval=x.^(-param.gamma);

for i=1:length(param.ygrid)
    
    cprime=funeval(c,fspace,[xprime,param.ygrid(i)*ones(ns,1)]);
     as=ns/length(param.ygrid);
      for j=1:length(param.ygrid)
        ypp(1+(j-1)*as:j*as,1)=param.yPP(j,i);
      end
    % ypp = param.ws(i);
    fval=fval-param.beta*(1+param.r)*ypp.*cprime.^(-param.gamma);
end

end

