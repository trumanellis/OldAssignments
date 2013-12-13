% Written by Truman Ellis
clear all; close all; clc;

nx = 100;
ny = 100;
x = linspace(-5, 3, nx);
y = linspace(-3, 5, ny);
U = zeros(nx, ny);

for i = 1:nx
    for j=1:ny
       U(i,j) = mueller(x(i), y(j)); 
    end
end

contour(x,y,exp(-.1*(U+200))',logspace(-9,2,100))%,linspace(-200,2000,500));
figure()
surf(x,y,exp(-.1*(U+200))')
colorbar
axis square
saveas(gcf,'MyMuller','pdf');