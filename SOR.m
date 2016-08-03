%------------SOR_New--------
clear all;
clc;
%-------constants------
M=0.5;
U=1;
ep=0.1;
E=10e-5;EP=1;
N=[100 150 200];
for g=1:3
w=1.25;
dx(g)=2*pi/(N(g));
dy(g)=2*pi/(N(g));
p=(dx(g)^2)/((1-M^2)*(dy(g)^2));
phi_old=zeros(N(g));
phi=zeros(N(g));
dphi_dx=zeros(N(g));
dphi_dy=zeros(N(g));
%-------initializing grid-----------
x=linspace(0,2*pi,N(g));
y=linspace(0,2*pi,N(g));
[X,Y] = meshgrid(x,y);
for i=1:N(g)
phi_old(i,1)=phi_old(i,2)-dy(g)*U*ep*cos((x(i)));
end
for i=1:N(g)
phi(i,1)=phi(i,2)-dy(g)*U*ep*cos((x(i)));
end
while EP>E
for i=2:N(g)-1
for j=2:N(g)-1
phi(i,j)=(1/(2*(1+p)))*(phi_old(i+1,j)+phi(i-1,j)
+p*(phi_old(i,j+1)+phi(i,j-1)));
end
end
phi(1,:)=phi(2,:);
phi(N(g),:)=phi(N(g)-1,:);
phi(:,N(g))=phi(:,N(g)-1);
for i=1:N(g)
phi(i,1)=phi(i,2)-dy(g)*U*ep*cos((x(i)));
endfor j=1:N(g)
for i=1:N(g)
phi(i,j)=(1-w)*phi_old(i,j)+w*phi(i,j);
%
phi_old(i,j)=phi(i,j);
end
end
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
str=sprintf('Contour plot of perturbation potential by SOR for grid size %i X
%i',N(g),N(g));
title(str);
for j=2:N(g)-1
for i=2:N(g)-1
dphi_dx(i,j) = ((phi(i+1,j)-phi(i-1,j))/2*dx(g));
dphi_dy(i,j)=((phi(i,j+1)-phi(i,j-1))/2*dy(g));
end
end
figure
hold on;
contourf(Y,X,dphi_dx,'edgecolor','none'); %contour plotxlabel('X');
ylabel('v');
title('Contour plot of v by SOR');
figure
hold on;
contourf(Y,X,dphi_dy,'edgecolor','none'); %contour plot
xlabel('X');
ylabel('u');
title('Contour plot of u by SOR');
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
