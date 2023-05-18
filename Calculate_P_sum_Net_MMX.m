clear all
clc
close all;
syms al be th r F real
assume((0<=th)&(th<=pi/2))
assume((0<=al)&(al<=1))
assume((0<=be)&(be<=1))
assume((0<=r)&(r<=1))
assume(al,'real')
assume(be,'real')
assume(th,'real')
assume(r,'real')
%assume(al^2+be^2==1)
%al=sqrt(1-be^2);
H = [1;0]; %|0>
V = [0;1]; %|1>
Id2 = eye(2,2);

%Maximally mixed initial state:
 si=al*kron(H,H)+be*kron(V,V);
 rnd=(1-F)*(1/4)*kron(Id2,Id2)+(F)*(si*si')
 rho_in=simplify(rnd/trace(rnd));

aqubit=[1,0;0,0];

CNOT=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0];
CNOT_R=[1,0,0,0;0,0,0,1;0,0,1,0;0,1,0,0];
M0=kron(H*H',kron(Id2,Id2));
M1=kron(V*V',kron(Id2,Id2));

k=8; %set k round

rho_combine=cell(1,k);
rho_success=cell(1,k);
rho_fail=cell(1,k);
p_success=cell(1,k);
p_total=cell(1,k);
 p_fail=cell(1,k);
 p_net=cell(1,k);
sum=0;
for i=1:k
   th=atan((al/be)^(2^(i-1))); 
   H=[cos(th),-sin(th);sin(th),cos(th)];
   if(i==1)
       rho_combine{i}=kron(CNOT_R,Id2)*kron(H*aqubit*H',rho_in)*kron(CNOT_R,Id2)';
   else
       rho_combine{i}=kron(CNOT_R,Id2)*kron(H*aqubit*H',rho_fail{i-1})*kron(CNOT_R,Id2)';
   end
   rho_success{i}=simplify(PartialTrace(M0*rho_combine{i}*M0',[1]));
   p_success{i}=simplify(trace(rho_success{i}));
   rf=simplify(PartialTrace(M1*rho_combine{i}*M1',[1]));
   rho_fail{i}=simplify(rf/trace(rf));
   p_fail{i}=trace(rf);
   mult=1;
   for jj=1:i-1
   mult=mult*p_fail{jj};
   end
   p_net{i}=p_success{i}*mult;
   sum=sum+p_net{i};
   p_total{i}=(sum);
end
celldisp( p_total)
%{
celldisp( p_success)
celldisp( p_fail)
celldisp(rho_success)
celldisp( p_fail)
celldisp( p_net)
%}
pp=1;
j=1;

for i=0.1:0.2:0.4
al=sqrt(i);%set initial state parameters
be=sqrt(1-al^2);

for k=0:0.3:1
    F=k;
 figure(pp);
hold on 
p_single_net_data=vpa(subs(p_success));
plot(p_single_net_data,'--o','LineWidth',2);

figure(pp+1);
hold on 
p_total_data=vpa(subs(p_total));
plot(p_total_data,'--o','LineWidth',2)
leg_str{j}=['F=',num2str(k)];
j=j+1;
end
figure(pp);
legend(leg_str)
title(['p\_success, \alpha^{2}=',num2str(i)])
figure(pp+1);
legend(leg_str)
title(['p\_sum\_net, \alpha^{2}=',num2str(i)])
pp=pp+2;
j=1;
end

