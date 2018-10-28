function F=mygpest(params)
% implement the gamma process based road crack modeling
%estimated parameters: 
%shape parameter v(t)=c t^b, b can be any fraction 1~2 with 0.1
%scale parameter (constant over time) mu
%params(1) for c, params(2) for mu

dattab1=table2array(readtable('trial1.txt','Delimiter','\t','ReadVariableNames',false));
datcurv=dattab1(4,1:end);
datpnum=length(datcurv);
b=1.1;
tval=linspace(0,datpnum-1,datpnum);
tpowvec=zeros(datpnum-1,1);
for i=2:datpnum
    tpowvec(i-1)=tval(i)^b-tval(i-1)^b;
end
F(1)=params(1)*tval(length(tval))^b/datcurv(end)-params(2);
psifunvec=zeros(datpnum-1,1);
for i=2:datpnum
   delta=datcurv(i)-datcurv(i-1);
   if delta==0
       delta=1e-6;
   end
   psifunvec(i-1)=psi(params(1)*tpowvec(i-1))-log(delta); 
end
F(2)=dot(tpowvec,psifunvec)-tval(length(tval))^b*log(params(1)*tval(length(tval))^b/datcurv(end));
end