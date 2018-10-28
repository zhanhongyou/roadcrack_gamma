clc
clear all

% single curves, MLE estimates
%paramest=fsolve(@mygpest,[1 5]);

% multiple curves, joint likelihood MLE estimates
b=linspace(1.1,2,10);
paramest_all=zeros(length(b),2);
%filnam='trial1.txt'; 
filnam='trial2.txt';
for k=1:length(b)
    myfunc = @(params) mygpest_multi(params,b(k),filnam);
    paramest_all(k,1:2)=fsolve(myfunc,[0.5 0.5]);
end

joinloglik=zeros(length(b),1);
dat1=importdata(filnam,'\t');
%calculate joint log-likelihood
for k=1:length(b)
    bval=b(k);
    cval=paramest_all(k,1);
    muval=paramest_all(k,2);
    for i=1:length(dat1(1:end,1))
        datcurv=dat1(i,~isnan(dat1(i,1:end)));
        datpnum=length(datcurv);

        tval=linspace(0,datpnum-1,datpnum);
        for j=2:datpnum
            xincval=datcurv(j)-datcurv(j-1);
            powval=cval*(tval(j)^bval-tval(j-1)^bval);
            joinloglik(k)=joinloglik(k)+powval*log(muval)-log(gamma(powval))+powval*log(xincval)-muval*xincval;
        end
    end
end
output1=horzcat(b',paramest_all,joinloglik);
[optval,optind]=max(output1(1:end,4));
optparam=output1(optind,1:end-1); %b,c,mu

%simulation using optimal parameters
% issues with simulation!!!!!!! shape function???
simcurvnum=10;
simtpnum=7;
tvec=linspace(0,simtpnum-1,simtpnum);
vtfunc=optparam(2)*(tvec.^optparam(1));%v(t)=ct^b
simxval=zeros(simcurvnum,simtpnum);
seedval=100;
for i=1:simcurvnum
   for j=2:simtpnum
      %rng(seedval)
      xinc=gamrnd((vtfunc(j)-vtfunc(j-1)),optparam(3));
      simxval(i,j)=simxval(i,j-1)+xinc;
      %seedval=seedval+1;
   end
   plot(tvec,simxval(i,1:end))
   hold on
end
hold off