clc
clear
close all

%poznate vrednosti 
A = 1; 
f0 = 0.2;
phi = 1;
sigma = 1;
N = 10;

% generisanje realizacija x
x = zeros(1,N);
for i=1:N
    x(1, i) = A*cos(2*pi*f0*(i-1) + phi) + randn*sigma;
end

%log verodostojnost i njen prikaz
phi_osa = -1:0.1:7;
suma = zeros(size(phi_osa));
for i=1:N
   suma = suma + (x(1,i) - A*cos(2*pi*f0*(i-1)+phi_osa)).^2;
end
l = -N/2*log(2*pi*sigma^2) - 1/(2*sigma^2)*suma;

figure(1)
plot(phi_osa,l);
title('Log-verodostojnost');
xlabel('\phi'); ylabel('l(\phi)');

%prvi izvod i njegov prikaz
suma = zeros(size(phi_osa));
for i=1:N
    suma = suma + x(1,i)*sin(2*pi*f0*(i-1)+phi_osa)-A/2*sin(4*pi*f0*(i-1)+2*phi_osa);
end
dl = -A/(2*sigma^2)*2*suma;
prvi = dl;
figure(2)
plot(phi_osa,dl);
title('Prvi izvod log-verodostojnosti');
xlabel('\phi'); ylabel('dl(\phi)/d\phi');

%drugi izvod i njegov prikaz
suma = zeros(size(phi_osa));
for i=1:N
    suma = suma + x(1,i)*cos(2*pi*f0*(i-1)+phi_osa)-A*cos(4*pi*f0*(i-1)+2*phi_osa);
end
ddl = -A/(2*sigma^2)*2*suma;

figure(3)
plot(phi_osa,ddl);
title('Drugi izvod log-verodostojnosti');
xlabel('\phi'); ylabel('d^2l(\phi)/d\phi^2');


%% drugi deo zadatka
pocetni = [0.5 1.5 4];
for i=1:3
    figure()
    hold on
    for k=1:5
        phi_k = pocetni(i);
        % generisanje realizacija x
        x = zeros(1,N);
        for j=1:N
            x(1, j) = A*cos(2*pi*f0*(j-1) + phi) + randn*sigma;
        end
        suma=0;
        for j=1:N
            suma = suma + x(1,j)*sin(2*pi*f0*(j-1)+phi_k)-A/2*sin(4*pi*f0*(j-1)+2*phi_k);
        end
        dl = -A/(2*sigma^2)*2*suma;
        suma = 0;
        for j=1:N
            suma = suma + x(1,j)*cos(2*pi*f0*(j-1)+phi_k)-A*cos(4*pi*f0*(j-1)+2*phi_k);
        end
        ddl = -A/(2*sigma^2)*2*suma;
        phi_err = phi_k - dl/ddl;
        phi_k1 = phi_err - 2*pi*floor(phi_err/2/pi);
        niz = [phi_k phi_k1];
        treshold = 0.01;
        br_iter = 9;
        while (abs(phi_k-phi_k1)>treshold && br_iter>0)
            br_iter = br_iter-1;
            phi_k=phi_k1;
            suma = 0;
            for j=1:N
                suma = suma + x(1,j)*sin(2*pi*f0*(j-1)+phi_k)-A/2*sin(4*pi*f0*(j-1)+2*phi_k);
            end
            dl = -A/(2*sigma^2)*2*suma;
            suma = 0;
            for j=1:N
                suma = suma + x(1,j)*cos(2*pi*f0*(j-1)+phi_k)-A*cos(4*pi*f0*(j-1)+2*phi_k);
            end
            ddl = -A/(2*sigma^2)*2*suma;
            phi_err = phi_k - dl/ddl;
            phi_k1 = phi_err - 2*pi*floor(phi_err/2/pi);
            niz(end+1) = phi_k1;
        end
        n = 0:length(niz)-1;
        plot(n,niz);
    end
    hold off
    xlabel('k');
    ylabel('phi_k');
    title(['NRM za phi_0 = ', num2str(pocetni(i))] );
    grid on
    legend('1. realizacija', '2. realizacija', '3. realizacija', '4. realizacija', '5. realizacija');
end

%% trazenje asimptotske FGV EMV
phi0 = 0.5;
sigma = sqrt(1);
N = 10;
Nr = 1000;
estimacija = zeros(1, Nr);
for j=1:Nr
    x = zeros(1, N);
    for i = 1:N
        x(1, i) = A*cos(2*pi*f0*(i-1) + phi) + randn*sigma;
    end
    phi_k1 = phi0;
    maxiter = 5;
    %kriterijum konvergencije - 5 iteracija
    while(maxiter > 0)
        maxiter = maxiter - 1;
        phi_k = phi_k1;
        suma = 0;
        for i=1:N
           suma = suma + x(1,i)*sin(2*pi*f0*(i-1) + phi_k) - A/2*sin(4*pi*f0*(i-1) + 2*phi_k); 
        end
        dl = -A/(sigma)^2*suma;
        suma = 0;
        for i=1:N
            suma = suma + x(1, i)*cos(2*pi*f0*(i-1) + phi_k) - A*cos(4*pi*f0*(i-1) + 2*phi_k); 
        end
        ddl = -A/(sigma)^2*suma;
        phierr = phi_k - dl/ddl;
        phi_k1 = phierr - 2*pi*floor(phierr/2/pi);
        %disp(phi_k1)
    end
    estimacija(1, j) = phi_k1;
end
figure();
histogram(estimacija, 25, 'Normalization', 'pdf');         
hold all
I = A^2*N/(sigma^2*2);
var = 1/I;
x = xlim;
x = x(1):0.01:x(2);
phi_a = 1/(sqrt(2*pi*var))*exp(-(1/(2*var)).*(x - phi).^2);
plot(x, phi_a);
legend('Estimirano','Teorijski');
hold off