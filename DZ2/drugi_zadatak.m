clc
clear
close all

FID = fopen('gps_data.txt');
formatSpec = '%f';
Opservacije = fscanf(FID, formatSpec);
%s=[p;v;a];
pr_put(1) = Opservacije(1);
for i=1:(length(Opservacije)-1)
pr_put(i+1) = Opservacije(i+1)-Opservacije(i);
end
ubrzanje(1) = pr_put(1);
for i=1:(length(pr_put)-1)
ubrzanje(i+1) = pr_put(i+1)-pr_put(i);
end
Ts = 1;
sgmw = 1;
sgmu = 0.002;
sgmu1 = 5*0.005;
sgmu2 = 0.005/5;
varu= sgmu^2;
varu1=sgmu1^2;
varu2=sgmu2^2;
H = [1, 0, 0];
B = [0; 0; 1];
A =[1, Ts, Ts^2/2; 0, 1, Ts; 0, 0, 1];
%Q = zeros(3,3);
Q = varu;
Q1=varu1;
Q2=varu2;
C = sgmw^2;
% -1 je 1
s_kapa = zeros(3,length(Opservacije));
s=[0;0;0];
%K_pojacanje(:,1)=[0;0;0];
%s_kapa(:,1) = [0;0;0];
s_kapa(:,1) = A*s;
M=eye(3);
M_pret(:,1) = diag(M);
M = A*M*transpose(A) + B*Q1*transpose(B);
M_tren(:,1) = diag(M);
K_pojacanje(:,1)=M*H'.*inv(C + H*M*H');
for i=2:(length(Opservacije)+1)
    
    s_kapa(:,i) = A*s_kapa(:,i-1);
    M = A*M*transpose(A) + B*Q1*transpose(B);
    M_pret(:,i) = diag(M);
    K_pojacanje(:,i) = M*H'.*inv(C + H*M*H');
    s_kapa(:,i) = s_kapa(:,i)+K_pojacanje(:,i)*(Opservacije(i-1)-H*s_kapa(:,i));
    M = (eye(3)-K_pojacanje(:,i)*H)*M;
    M_tren(:,i) = diag(M);
end
figure(1)
plot(1:length(s_kapa),s_kapa(3,:));
title('Procenjeno ubrzanje');
figure(2)
plot(1:length(s_kapa),s_kapa(2,:));
title('Procenjena brzina');

figure(3)
plot(1:length(s_kapa),s_kapa(1,:));
hold all
plot(1:length(Opservacije),Opservacije);
title('Pozicije');
legend('estimirano','opservirano');
hold off


figure(4)
plot(1:length(K_pojacanje),K_pojacanje);
title('Pojacanje K');
legend('pozicija','brzina','ubrzanje');

figure(5)
plot(1:length(M_pret),M_pret);
title('M[n|n-1]');
legend('pozicija','brzina','ubrzanje');

figure(6)
plot(1:length(M_tren),M_tren);
title('M[n|n]');
legend('pozicija','brzina','ubrzanje');

%% za 5sigma

M = A*M*transpose(A) + B*Q1*transpose(B);
%M_tren(:,1)=[1;1;1];
M_tren(:,1) = diag(M);
K_pojacanje(:,1)=M*H'.*inv(C + H*M*H');
for i=2:(length(Opservacije))
    
    s_kapa(:,i) = A*s_kapa(:,i-1);
    M = A*M*transpose(A) + B*Q1*transpose(B);
    M_pret(:,i) = diag(M);
    K_pojacanje(:,i) = M*H'.*inv(C + H*M*H');
    s_kapa(:,i) = s_kapa(:,i)+K_pojacanje(:,i)*(Opservacije(i)-H*s_kapa(:,i));
    M = (eye(3)-K_pojacanje(:,i)*H)*M;
    M_tren(:,i) = diag(M);
end
figure()
plot(1:length(s_kapa),s_kapa(3,:));
title('Procenjeno ubrzanje,\sigma=5\sigma_{opt}');
figure()
plot(1:length(s_kapa),s_kapa(2,:));
title('Procenjena brzina,\sigma=5\sigma_{opt}');

figure()
plot(1:length(s_kapa),s_kapa(1,:));
hold all
plot(1:length(Opservacije),Opservacije);
title('Pozicije,\sigma=5\sigma_{opt}');
legend('estimirano','opservirano');
hold off


figure()
plot(1:length(K_pojacanje),K_pojacanje);
title('Pojacanje K');
legend('pozicija','brzina','ubrzanje');

figure()
plot(1:length(M_pret),M_pret);
title('M[n|n-1]');
legend('pozicija','brzina','ubrzanje');

figure()
plot(1:length(M_tren),M_tren);
title('M[n|n]');
legend('pozicija','brzina','ubrzanje');

%% sigma=sigma/5

M = A*M*transpose(A) + B*Q2*transpose(B);
%M_tren(:,1)=[1;1;1];
M_tren(:,1) = diag(M);
K_pojacanje(:,1)=M*H'.*inv(C + H*M*H');
for i=2:(length(Opservacije))
    
    s_kapa(:,i) = A*s_kapa(:,i-1);
    M = A*M*transpose(A) + B*Q2*transpose(B);
    M_pret(:,i) = diag(M);
    K_pojacanje(:,i) = M*H'.*inv(C + H*M*H');
    s_kapa(:,i) = s_kapa(:,i)+K_pojacanje(:,i)*(Opservacije(i)-H*s_kapa(:,i));
    M = (eye(3)-K_pojacanje(:,i)*H)*M;
    M_tren(:,i) = diag(M);
end
figure()
plot(1:length(s_kapa),s_kapa(3,:));
title('Procenjeno ubrzanje,\sigma=\sigma_{opt}/5');
figure()
plot(1:length(s_kapa),s_kapa(2,:));
title('Procenjena brzina,\sigma=\sigma_{opt}/5');

figure()
plot(1:length(s_kapa),s_kapa(1,:));
hold all
plot(1:length(Opservacije),Opservacije);
title('Pozicije,\sigma=sigma_{opt}/5');
legend('estimirano','opservirano');
hold off


figure(4)
plot(1:length(K_pojacanje),K_pojacanje);
title('Pojacanje K');
legend('pozicija','brzina','ubrzanje');

%% promenjeni pocetni uslovi
s_pocetno = [5;5;5];
M_pocetno = 10*eye(3);
s = [5;5;5];
M=10*eye(3);



M_pret(:,1) = diag(M);
M = A*M*transpose(A) + B*Q*transpose(B);
M_tren(:,1) = diag(M);
K_pojacanje(:,1)=M*H'.*inv(C + H*M*H');
for i=2:(length(Opservacije))
    
    s_kapa(:,i) = A*s_kapa(:,i-1);
    M = A*M*transpose(A) + B*Q*transpose(B);
    M_pret(:,i) = diag(M);
    K_pojacanje(:,i) = M*H'.*inv(C + H*M*H');
    s_kapa(:,i) = s_kapa(:,i)+K_pojacanje(:,i)*(Opservacije(i)-H*s_kapa(:,i));
    M = (eye(3)-K_pojacanje(:,i)*H)*M;
    M_tren(:,i) = diag(M);
end
figure()
plot(1:length(s_kapa),s_kapa(3,:));
title('Procenjeno ubrzanje');
figure()
plot(1:length(s_kapa),s_kapa(2,:));
title('Procenjena brzina');

figure()
plot(1:length(s_kapa),s_kapa(1,:));
hold all
plot(1:length(Opservacije),Opservacije);
title('Pozicije');
legend('estimirano','opservirano');
hold off


figure()
plot(1:length(K_pojacanje),K_pojacanje);
title('Pojacanje K');
legend('pozicija','brzina','ubrzanje');

figure()
plot(1:length(M_pret),M_pret);
title('M[n|n-1]');
legend('pozicija','brzina','ubrzanje');

figure()
plot(1:length(M_tren),M_tren);
title('M[n|n]');
legend('pozicija','brzina','ubrzanje');
