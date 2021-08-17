clc
clear
close all

N = 10^5;

%definisanje parametara
m1 = 0; s1 = 1;
m2 = 0; s2 = 0.5;

rho = 0;
m = [m1 m2];
sigma = [s1^2 rho*s1*s2; rho*s1*s2 s2^2];


%generisanje N odbiraka slucajnog vektora X i eksperimentalna fgv
u = randn(N,2);
x_eksp(:,1) = m(1) + u(:,1).*(s1);
x_eksp(:,2) = m(2) + u(:,2).*(s2);
figure(1)
hst = histogram2(x_eksp(:,1),x_eksp(:,2),20,'Normalization','pdf');
xlabel('x_1');ylabel('x_2');zlabel('$$\hat{f}(x_1,x_2)$$','Interpreter','Latex')
title('Eksperimentalno procenjena fgv');

%analiticki odredjena fgv
x1 = -3*s1:0.1:3*s1;
x2 = -3*s2:0.1:3*s2;

[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];
F = 1/(2*pi*s1*s2)*exp(-(X1.^2)/2/(s1^2)-(X2.^2)/2/(s2^2));

figure(2)
surf(x1,x2,F);
xlabel('x_1'); ylabel('x_2'); zlabel('f(x_1,x_2)');
title('Analiticki odredjena fgv');

%generisanje odbiraka za Y
A = [1.6, -2.4; 0, 6];
b = [0; 10];
Y = A*x_eksp' + b;
Y = Y';

%vektor ocekivanja
m = sum(Y)/N;
disp(['Vektor ocekivanja EY']);
disp(m);

%kovarijaciona matrica
% varY1 = sum((Y(:,1)-m(1,1)).^2)/(N-1);
% varY2 = sum((Y(:,2)-m(1,2)).^2)/(N-1);
% cov_Y1Y2 = sum((Y(:,1)-m(1,1)).*(Y(:,2)-m(1,2)))/(N-1);
% R = [varY1 cov_Y1Y2; cov_Y1Y2 varY2];

R = ((Y-m)'*(Y-m))/(N-1);
disp(['Kovarijaciona matrica:']);
disp(R);

%ro = 0
Y1 = [2*x_eksp(:, 1), 6*x_eksp(:, 2)+10];

%ro = 0.8
Y2 = [1.2*x_eksp(:, 1)+3.2*x_eksp(:, 2), 6*x_eksp(:, 2)+10];

%crtanje grafika
figure(3);
plot(Y(:, 1), Y(:, 2), 'x');
xlabel('Y_1');
ylabel('Y_2');
title('Odbirci vektora Y kada je \rho = -0.6');

figure(4);
plot(Y1(:, 1), Y1(:, 2), 'x');
xlabel('Y_1');
ylabel('Y_2');
title('Odbirci vektora Y kada je \rho = 0');


figure(5);
plot(Y2(:, 1), Y2(:, 2), 'x');
xlabel('Y_1');
ylabel('Y_2');
title('Odbirci vektora Y kada je \rho = 0.8');





