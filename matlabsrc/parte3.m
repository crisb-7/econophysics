clc; clear;

N = 1000;    M = 1e+5;       C = 50;    alpha = 1/300;
mmax = (M/(sum(1:C)*(N/C))*C)*2; 
Class = linspace(0, 650, C+1);
Mk = linspace(0, mmax, C);
Clases = zeros(5, N); % Matriz que guarda los estados de las clases
N_vec = 1:N;        V = 0;      Y = 12;     wa = 10;     wb = 90;
avgW = (wa + wb)/2;
myCenters = 0.5*( Class(1:C) + Class(2:C+1) );

%% Simulacion 1: Distribucion delta

Clases(3,:) = ones(1,N);
Clases(4,:) = ones(1,N)*(M/N); %Distribucion Delta
ny = 1;   X = zeros(1,ny);       W = zeros(1,ny);
%T = 12;
T = 12*100;         % T = ny*Y;
t = linspace(0,T, T);
valEsp = zeros(1,C);
for i=1:T
    %{
    if mod(i,50) == 0
        histogram(Clases(4,Clases(1,:)~=0), Class)
        grid on
        hold on
        histogram(Clases(4,Clases(2,:)~=0), Class)
        histogram(Clases(4,Clases(3,:)~=0), Class)
        set(gca, "FontSize", 13, "FontName", "Times")
        xlabel("Dinero (\$)", "Interpreter", "latex")
        ylabel("N\'umero de personas ($n_k$)", "Interpreter", "latex")
        text = sprintf("Distribuci\'on del dinero, t = %d mes(es)", i);
        legend("Capitalistas", "Trabajadores", "Desempleados", "Interpreter", "latex")
        %title("Distribuci\'on del dinero", "Interpreter", "latex")
        title(text, "Interpreter", "latex")
        pause(1/120)
        hold off
    end
    %}
    x = 0;
    dobleU = 0;
    for j=1:N
        a = randi(N);
        Clases = HiringRule(a, Clases, avgW);
        b = randi(N);
        [Clases, V] = ExpenditureRule(b, a, N, Clases, V);
        [Clases, V, x, X] = MarketSample(a, Clases, V, i, x, X);
        Clases = FiringRule(Clases, a, avgW);
        [Clases, dobleU, W] = WagePay(a, Clases, wa, wb, i, dobleU, W);
    end
    if mod(i,12) == 0
        incT = sort(Clases(5,:)); 
        incC = sort(Clases(5, Clases(1,:)~=0));
        incW = Clases(5,Clases(2,:)~=0);
        %Clases(5,:) = zeros(1,N);   
    end
    MAX(i) = max(Clases(4,:));
    nk = histcounts(Clases(4,:), Class);
    S(i) = N*log(N) - sum(nk.*log(nk), "omitnan");
    OcA(i) = sum(alpha*myCenters.*nk);
    OcB(i) = sum(nk.*(1-exp(-alpha*myCenters)));
    
    VEm = 0;
    for u = 1:C
        VEm = VEm + myCenters(u).*nk(u)/N;
    end
    VE(i) = VEm;
end

n = 10*N;
CBin1 = linspace(0, max(MAX), n); % Ingresos
CBin2 = linspace(0, M, n); % Dinero circulante

%% Graficas Dinero Circulante

% --------------------- Distribucion completa --------------------------
%{
%mT = Clases(4,Clases(4,:) > 0);
mT = Clases(4,:);
hmT = histogram(mT,CBin1, "Normalization", "pdf");
cdfmT = hmT.Values;
cdfmT_vec = 0.5*(hmT.BinEdges(1:n-1) + hmT.BinEdges(2:n)); 
%cdfmT = cdfmT(1-cdfmT > 0.00001);
%cdfmT_vec = cdfmT_vec(1-cdfmT > 0.00001);
semilogx(cdfmT_vec, cdfmT, "b*-", "LineWidth", 1.2)
grid on
xlabel("Dinero circulante (\$)","Interpreter", "latex")
ylabel("Probabilidad de ingresos $P(M>m)$", "Interpreter", "latex")
title("Distribucion Completa","Interpreter", "latex")
set(gca, "FontSize", 13, "FontName", "Times")
legend("Poblaci\'on total","Interpreter", "latex")
%}

%---------------- Distribucion Desagregada ------------------
%{
mC = Clases(4,Clases(1,:) ~= 0);        mC = mC(mC > 0 );
mW = Clases(4,Clases(2,:) ~= 0);        mW = mW(mW > 0 );
hC = histogram(mC,CBin1, "Normalization", "cdf");
cdfC = hC.Values;
cdfC_vec = 0.5*(hC.BinEdges(1:n-1) + hC.BinEdges(2:n)); 
cdfC = cdfC(1-cdfC > 0.0000001);
cdfC_vec = cdfC_vec(1-cdfC > 0.0000001);
hW = histogram(mW,CBin1, "Normalization", "cdf");
cdfW = hW.Values;
cdfW_vec = 0.5*(hW.BinEdges(1:n-1) + hW.BinEdges(2:n)); 
cdfW = cdfW(1-cdfW > 0.0000001);
cdfW_vec = cdfW_vec(1-cdfW > 0.0000001);
loglog(cdfC_vec, 1-cdfC, "b*-", "LineWidth", 1)
grid on
hold on
loglog(cdfW_vec, 1-cdfW, "r*-", "LineWidth", 1)
xlabel("Dinero circulante (\$)","Interpreter", "latex")
ylabel("Probabilidad de dinero $P(M>m)$", "Interpreter", "latex")
title("Dinero Circulante Desagregado","Interpreter", "latex")
set(gca, "FontSize", 13, "FontName", "Times")
legend("Capitalistas", "Trabajadores", "Interpreter", "latex")
%}

%% Graficas Ingresos

% --------------------- Distribucion completa --------------------------
%incomeT = incT(incT > 0);
incomeC = incC(incC > 0);
incomeW = incW(incW > 0);

% incomeT = Clases(5,:);
% incomeC = Clases(5, Clases(1,:)~=0);
% incomeW = Clases(5, Clases(2,:)~=0);

%{
hT = histogram(incT,CBin1, "Normalization", "cdf");
cdfT = hT.Values;
cdfT_vec = 0.5*(hT.BinEdges(1:n-1) + hT.BinEdges(2:n)); 
% cdfT = cdfT(1-cdfT > 0.0001);
% cdfT_vec = cdfT_vec(1-cdfT > 0.0001);
loglog(cdfT_vec, 1-cdfT, "b*-", "LineWidth", 1.2)
grid on
xlabel("Ingresos (\$)","Interpreter", "latex")
ylabel("Probabilidad de ingresos $P(M>m)$", "Interpreter", "latex")
title("Distribucion Completa","Interpreter", "latex")
set(gca, "FontSize", 13, "FontName", "Times")
legend("Poblaci\'on total","Interpreter", "latex")
%}

% ---------------- Distribucion Desagregada ------------------

%{
hC = histogram(incomeC,CBin1, "Normalization", "cdf");
cdfC = hC.Values;
cdfC_vec = 0.5*(hC.BinEdges(1:n-1) + hC.BinEdges(2:n)); 
cdfC = cdfC(1-cdfC > 0.0001);
cdfC_vec = cdfC_vec(1-cdfC > 0.0001);
hW = histogram(incomeW,CBin1, "Normalization", "cdf");
cdfW = hW.Values;
cdfW_vec = 0.5*(hW.BinEdges(1:n-1) + hW.BinEdges(2:n)); 
cdfW = cdfW(1-cdfW > 0.0001);
cdfW_vec = cdfW_vec(1-cdfW > 0.000001);
loglog(cdfC_vec, 1-cdfC, "b*-")
grid on
hold on
loglog(cdfW_vec, 1-cdfW, "r*-")
xlabel("Ingresos (\$)","Interpreter", "latex")
ylabel("Probabilidad de ingresos $P(M>m)$", "Interpreter", "latex")
title("Distribucion Desagregada","Interpreter", "latex")
set(gca, "FontSize", 13, "FontName", "Times")
legend("Capitalistas", "Trabajadores", "Interpreter", "latex")
%}

%% Graficas Bienestar

%
figure()
plot(t,OcA, "r", "LineWidth", 1.3)
hold on
plot(t,OcB, "b", "Linewidth", 1.3)
grid on
axis square
set(gca, "FontSize", 13, "FontName", "Times")
xlabel("Tiempo $t$", "Interpreter", "latex")
ylabel("Bienestar Colectivo $O(n_k)$", "Interpreter", "latex")
title("Funcion Objetivo", "Interpreter", "latex")
legend("A) Wright", "B) Wright", "Interpreter", "latex")
axis square
%}

%% Graficas Entropia

%{
figure()
plot(t, S, "b", "LineWidth", 1.3)
grid on
axis square
set(gca, "FontSize", 13, "FontName", "Times")
xlabel("Tiempo $t$", "Interpreter", "latex")
ylabel("Entrop\'ia $S(t)$ (meses)", "Interpreter", "latex")
title("Entrop\'ia en funci\'on del tiempo", "Interpreter", "latex")
axis square
%}

function Clases = HiringRule(a, Clases, avgW)
    if Clases(2,a) == 0 && Clases(1,a) == 0 %Si agente a es desempleado
        P = zeros(1,length(Clases));
        smC = sum(Clases(4, Clases(1,:) == 1));
        P(Clases(1,:)==1) = Clases(4, Clases(1,:) == 1)/smC;
        P(Clases(3,:)==1) = Clases(4, Clases(3,:) == 1)/smC;
%         plot(P)
%         pause(1/60)
        while Clases(2,a) == 0  % Buscar empleador hasta encontrar
            c = randsample(length(Clases),1, true, P); % Posible fuente de error
            if Clases(4,c) > avgW && c ~= a
                Clases(1,a) = 0;
                Clases(2,a) = c;        % Cambiar estado a trabajador
                Clases(3,a) = 0;
                
                Clases(1,c) = 1;
                Clases(2,c) = 0;
                Clases(3,c) = 0;
            end
        end
        %P = zeros;
    end
end

function [Clases, V] = ExpenditureRule(b, a, N, Clases, V)
    while b == a
        b = randi(N);
    end
    m = (Clases(4,b)*rand());
    %m = randi([0, Clases(4,b)],1);
    V = V + m;
    Clases(4,b) = Clases(4,b) - m;
end

function [Clases, V, x, X] = MarketSample(a, Clases, V, i, x, X)
    %m = randi([0 V]);
    m = rand()*V;
    V = V - m;
    if Clases(2,a) ~= 0 && Clases(1,a) == 0
        ea = Clases(2,a);
        Clases(4, ea) = Clases(4, ea) + m;
        Clases(5, ea) = Clases(5, ea) + m;
    end
    if Clases(1,a) ~= 0 && Clases(2,a) == 0
        Clases(4,a) = Clases(4,a) + m;
        Clases(5,a) = Clases(5,a) + m;
    end
    x = x + m;
    if mod(i,12) == 0
        X(i/12) = x;
    end
end

function Clases = FiringRule(Clases, a, avgW)
    if Clases(1,a) == 1
        ma = Clases(4,a);
        F = find(Clases(2,:)==a);
        nEmp = sum(Clases(2, :) == a);
        u = max( nEmp - floor(ma./avgW) ,0);
        if u > 0
            f = randsample(F,u);
            Clases(1,f) = 0;
            Clases(2,f) = 0;
            Clases(3,f) = 1;
            if length(F) == u
                Clases(1, a) = 0;
                Clases(3, a) = 1;
            end
        end
    end
end

function [Clases, dobleU, W] = WagePay(a, Clases, wa, wb, i, dobleU, W)
    if Clases(1,a) == 1
        Wa = find(Clases(2,:) == a);
        for j=Wa
            if Clases(4,a) > wb
                w = (wb-wa)*rand + wa;
            elseif Clases(4,a) > 0
                w = rand()*Clases(4,a);
            else
                w = 0;
            end
            Clases(4,j) = Clases(4,j) + w;
            Clases(5,j) = Clases(5,j) + w;
            Clases(4,a) = Clases(4,a) - w;
            
            dobleU = dobleU + w;
            if mod(i,12) == 0
                W(i/12) = dobleU;
            end
        end
        if Clases(4,a) < 0
            disp("WagePay")
        end
    end
end