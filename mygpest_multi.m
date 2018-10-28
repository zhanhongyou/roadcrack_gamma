function F=mygpest_multi(params,b,filnam)
% implement the gamma process based road crack modeling
%estimated parameters: 
%shape parameter v(t)=c t^b, b can be any fraction 1~2 with 0.1
%scale parameter (constant over time) mu
%params(1) for c, params(2) for mu

dattab1=importdata(filnam,'\t');

item1=0;
item3=0;
tpowtot=0;
xdifftot=0;
for i=1:length(dattab1(1:end,1))
    datcurv=dattab1(i,~isnan(dattab1(i,1:end)));
    datpnum=length(datcurv);

    tval=linspace(0,datpnum-1,datpnum);
    for j=2:datpnum
        tpowval=tval(j)^b-tval(j-1)^b;
        delta=datcurv(j)-datcurv(j-1);
        psival=psi(params(1)*tpowval);
        item1=item1+tpowval*(psival-log(params(1)));
        tpowtot=tpowtot+tpowval;
        xdifftot=xdifftot+delta;
        item3=item3+tpowval*log(delta);
    end 
end

F(1)=params(1)*tpowtot/xdifftot-params(2);
F(2)=item1-tpowtot*log(tpowtot/xdifftot)-item3;
end