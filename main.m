t=0;
Tin1=298;Tm=298;
data=[];
while t<180.1
q=1000;Cp=4180;qm1=0.05;hs=10000;
As=0.01;h1=5000;A1=0.1;M=0.3;
Tout1=q/(Cp*qm1)+Tin1;
Ts=q/(hs*As)+(Tin1+Tout1)/2;
Tin1=(Cp*qm1-h1*A1/2)/(h1*A1/2+Cp*qm1)*Tout1+h1*A1/(h1*A1/2+Cp*qm1)*Tm;
Cpm=fCpm(Tm);
data=[data;t,Tm,Tin1,Tout1,Ts];
deltat=0.1;
dTmdt = Cp*qm1*(Tout1-Tin1)/(Cpm*M);
Tm=Tm+dTmdt*deltat;
t=deltat+t;
end
figure;
plot(data(:,1),data(:,2:end));
legend('Tm','Tin1','Tout1','Ts')
xlabel('t');
