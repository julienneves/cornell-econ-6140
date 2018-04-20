function fval = euler_ca(x,c,fspace,s,param)

fval=x.^(-param.gamma);
as=size(s,1)/length(param.ygrid);
    
for i=1:length(param.ygrid)    
    cprime=funeval(c,fspace,(1+param.r)*(s-x)+exp(param.ygrid(i)));
%     for j=1:length(param.ygrid)
%         ypp(1+(j-1)*as:j*as,1)=param.yPP(j,i);
%     end
    ypp = param.ws(i);
    
    fval=fval-param.beta*(1+param.r)*ypp.*cprime.^(-param.gamma);
end

end

