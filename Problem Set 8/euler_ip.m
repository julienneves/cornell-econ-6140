function fval = euler_ip(x,c,fspace,s,param)

ns=size(s,1);
a=s(:,1);
y=exp(s(:,2));
aprime=(1+param.r)*a+y-x;
fval=x.^(-param.gamma);
as=ns/length(param.ygrid);

for i=1:length(param.ygrid)
  cprime=funeval(c,fspace,[aprime,param.ygrid(i)*ones(ns,1)]);
  for j=1:length(param.ygrid)
    ypp(1+(j-1)*as:j*as,1)=param.yPP(j,i);
  end
  fval=fval-param.beta*(1+param.r)*ypp.*cprime.^(-param.gamma);
end
