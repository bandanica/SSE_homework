clc
clear
close all


N = 100000;
broj = 50;

%generisanje slucajne promenljive Y
x = rand(N,1);
Y = zeros(N,1);
Y(x>=0 & x<0.5) = sqrt(2*x(x>=0 & x<0.5)) - 2;
Y(x>=0.5 & x<1) = -sqrt(2-2*x(x>=0.5 & x<1));

%prikaz fgv
figure(1)
hst = histogram(Y,broj, 'Normalization','pdf');
hold all
title('funkcija gustine verovatnoce'); xlabel('y'); ylabel('f_Y(y)');
%analiticka fgv-prikaz
y = [-2 -1 0 1];
fy = [0 1 0 0];
plot(y,fy, 'LineWidth',2);
xlim([-2.1 0.1]);
legend('procena f(y)','f(y)');

%racunanje varijanse i ocekivanja
ocekivanje_analiticki = -1;
varijansa_analiticki = 1/6;
ocekivanje_eksp = sum(Y)/N;
varijansa_eksp = sum((Y-ocekivanje_eksp).^2)/(N-1);
