% MAE 542 Midterm Project
%(1-M^2)phi_xx+phi_yy=0
%
%Analytical Solution
clear all;
clc;
%-------constants------
M=0.5;
U=1;
L=2*pi;
ep=0.1;
N=100;
dx=L/2*pi;
dy=L/2*pi;
dphi_dx=zeros(N);dphi_dy=zeros(N);
%-------initializing grid-----------
x=linspace(0,2*pi);
y=linspace(0,2*pi);
[X,Y] = meshgrid(x,y);
%-------------analytical solution-------------
phi=((-U.*ep).*exp(-(sqrt(1-M^2).*Y)).*cos(X))./sqrt(1-M^2);
figure
hold on;
contourf(X,Y,phi,'edgecolor','none'); %contour plot
for j=2:N-1
for i=2:N-1
dphi_dx(i,j) = ((phi(i+1,j)-phi(i-1,j))/2*dx);
dphi_dy(i,j)=((phi(i,j+1)-phi(i,j-1))/2*dy);
end
end
figure
hold on;
contourf(X,Y,dphi_dx,'edgecolor','none'); %contour plot
figure
hold on;
contourf(X,Y,dphi_dy,'edgecolor','none'); %contour plot
