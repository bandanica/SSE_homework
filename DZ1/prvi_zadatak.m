clc
clear 
close all

br_realizacija = 1000;

F = [0 0.01 0.34 1];
p = [0.01 0.33 0.66];
t = [0 1 2];

%generisanje slucajne promenljive
u = rand(br_realizacija,1);
x = zeros(br_realizacija,1);
x(u<=p(1)) = 0;
x(u>p(1) & u<=(p(1)+p(2))) = 1;
x(u>(p(1)+p(2))) = 2;

%prikaz histograma
figure(1)
C = categorical(x,[0 1 2],{'ivica','pismo','glava'});
h = histogram(C,'BarWidth',0.5);
title('Histogram za 10^3 ishoda'); ylabel('broj ishoda');

%procena funkcije mase verovatnoce
p_kapica(1) = numel(x(x==0))/br_realizacija;
p_kapica(2) = numel(x(x==1))/br_realizacija;
p_kapica(3) = numel(x(x==2))/br_realizacija;

%racunanje funkcije raspodele
F_kapica(1) = 0;
for i=1:length(p_kapica)
    F_kapica(i+1) = sum(p_kapica(1:i));
end

%prikaz analiticki i eksperimentalno dobijene
%funkcije raspodele
figure(2)
stairs([-1 0 1 2 3],[F, F(4)]);
hold all
title('Funkcije raspodele Fx(k)'); xlabel('k');ylabel('F_x(k)');
stairs([-1 0 1 2 3],[F_kapica, F_kapica(4)]);
leg = legend('$F_x(k)$','$\hat{F}_x(k)$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',12);
xlim([-0.5 2.5]);
ylim([-0.2 1.5]);
grid on

%prikaz analiticki i eksperimentalno dobijene
%funkcije mase verovatnoce
figure(3)
stem(0:length(p)-1, p);
hold all
stem(0:length(p_kapica)-1, p_kapica);
title('Funkcija mase verovatnoce');
xlabel('k');
ylabel('p_x(k)');
leg = legend('$p_x(k)$','$\hat{p}_x(k)$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',12);
axis([-0.5 2.5 0 1]);

%odredjivanje EX i var
varijansa = zeros(2,1);
ocekivanje = zeros(2,1);
ocekivanje(1) = 1.65;
varijansa(1) = 0.2475;


ocekivanje(2) = sum(x)/br_realizacija;
varijansa(2) = sum((x-ocekivanje(2)).^2)/(br_realizacija-1);

tabela = table(ocekivanje, varijansa);
disp(tabela);

