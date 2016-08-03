% MAE 542 Midterm Project
%(1-M^2)phi_xx+phi_yy=0
%
%Solution Using ADI
clear all;
clc;
%-------constants------
M=0.5;
U=1;
L=2*pi;ep=0.1;
EX=10e-5;
EP=1;
N=[100 150 200];
for g=1:3
dx(g)=2*pi/(N(g));
dy(g)=2*pi/(N(g));
p=(dx(g)^2)/((1-M^2)*(dy(g)^2));
phi_old=zeros(N(g));
phi=zeros(N(g));
phi_new=zeros(N(g));
%-------initializing grid-----------
x=linspace(0,2*pi,N(g));
y=linspace(0,2*pi,N(g));
[X,Y] = meshgrid(x,y);
while EP>EX
for j=1:N(g)
A(1,j)=0;
B(1,j)=-1;
C(1,j)=1;
D(1,j)=0;
end
for i=2:N(g)-1
for j=2:N(g)-1
A(i,j)=0.75/(dx(g)*dx(g));
B(i,j)=-2*((0.75/(dx(g)*dx(g)))+1/(dy(g)*dy(g)));
C(i,j)=0.75*(1/(dx(g)*dx(g)));
D(i,j)=-(dy(g)^(-2))*(phi_old(i,j+1)+phi(i,j-1));
end
end
for j=1:N(g)
A(N(g),j)=-1;
B(N(g),j)=1;C(N(g),j)=0;
D(N(g),j)=0;
end
TRI_2D_X(1,N(g),2,N(g)-1,A,B,C,D);
for j=1:N(g)
for i=1:N(g)
phi(i,j)=D(i,j);
phi_old(i,j)=phi(i,j);
end
end
phi(1,1)=phi(2,1);
phi(N(g),1)=phi(N(g)-1,1);
phi(1,N(g))=phi(2,N(g));
phi(N(g),N(g))=phi(N(g)-1,N(g));
% -----------------------------------------------
for i=1:N(g)
E(i,1)=0;
F(i,1)=-1;
G(i,1)=1;
H(i,1)=2*dy(g)*0.1*cos(x(i));
end
for i=2:N(g)-1
for j=2:N(g)-1
E(i,j)=(dy(g)^(-2));
F(i,j)=-2*((0.75*(1/(dx(g)*dx(g))))+1/(dy(g)*dy(g)));
G(i,j)=(1/(dy(g)*dy(g)));
H(i,j)=-(0.75*(1/(dx(g)*dx(g))))*(phi(i,j+1)+phi_new(i,j-1));
end
end
for i=1:N(g)
E(i,N(g))=-1;
F(i,N(g))=1;G(i,N(g))=0;
H(i,N(g))=0;
end
TRI_2D_Y(2,N(g)-1,1,N(g),E,F,G,H);
for j=1:N(g)
for i=1:N(g)
phi_new(i,j)=H(i,j);
phi(i,j)=phi_new(i,j);
end
end
phi_new(1,1)=phi_new(2,1);
phi_new(N(g),1)=phi_new(N(g)-1,1);
phi_new(1,N(g))=phi_new(2,N(g));
phi_new(N(g),N(g))=phi_new(N(g)-1,N(g));
for j=1:N(g)
for i=1:N(g)
phi_dummy(i,j)=phi_old(i,j);
phi_old(i,j)=phi(i,j);
end
end
s=abs(phi-phi_dummy);
EP=(sum(s(:)))/sum(abs(phi_dummy(:)));
end
%-------plot-----------------------
figure
hold on;
contourf(Y,X,phi,'edgecolor','none'); %contour plot
xlabel('X');
ylabel('Phi');
str=sprintf('Contour plot of perturbation potential by ADI for grid size %i X
%i',N(g),N(g));title(str);
for j=2:N(g)-1
for i=2:N(g)-1
dphi_dx(i,j) = ((phi(i+1,j)-phi(i-1,j))/2*dx(g));
dphi_dy(i,j)=((phi(i,j+1)-phi(i,j-1))/2*dy(g));
end
end
figure
hold on;
contourf(Y,X,dphi_dx,'edgecolor','none'); %contour plot
xlabel('X');
ylabel('v');
title('Contour plot of v by ADI');
figure
hold on;
contourf(Y,X,dphi_dy,'edgecolor','none'); %contour plot
xlabel('X');
ylabel('u');
title('Contour plot of u by ADI');
phi_comp=((-U.*ep).*exp(-(sqrt(1-M^2).*Y)).*cos(X))./sqrt(1-M^2);
for i=1:N(g)
for j=1:N(g)
rms(i,j)=(phi(i,j)-phi_comp(i,j))^2;
end
end
r=sum(rms(:));
error(g)=(1/(N(g)*N(g)))*sqrt(r);
end
figure;
hold on;
plot(log(dx),log(error));
xlabel('log(dx)');
ylabel('log(error)');
title('Log-Log plot of Error Vs dx');
