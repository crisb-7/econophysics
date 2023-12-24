%% Parámetros iniciales

clc; clear;
N = 1000;    M = 1e+5;        C = 50;   mmax = (M/(sum(1:C)*(N/C))*C)*2;
%mmax = 650;
Class = linspace(0, mmax, C+1); 
T = 12000;        Ci = length(Class);     dm = M*4e-5;
lambda = 0;      tau = 30;      a = 1/300;
Mk = linspace(0, mmax, C);

%% Caso 1: Distribución Delta Inicial

% t = 0, 50, 100, 200, 300
myCenters = 0.5*( Class(1:C) + Class(2:C+1) );
Mk_vec = linspace(0, mmax, N);
iva = 0;
mD = ones(1,N)*(M/N); % Distribucion Delta
Sdelta = zeros(1,T);
for t=1:T
	%{
    histogram(mD, Class, "FaceColor", "b")
    ylim([0, 150])
    grid on
    set(gca, "FontSize", 13, "FontName", "Times")
    xlabel("Dinero (\$)", "Interpreter", "latex")
    ylabel("N\'umero de personas ($n_k$)", "Interpreter", "latex")
    %text = sprintf("Distribuci\'on del dinero, t = %d", t);
    %legend(text, "Interpreter", "latex")
    title("Distribuci\'on del dinero", "Interpreter", "latex")
    pause(1/120)
    %}
    
    if mod(t,tau) == 0
        mD = mD + iva/N;
        iva = 0;
    end
    for j=1:N
        l = randi([1 N],1,1);
        s = sign(2*rand()-1);
        if mD(j) + dm*s < 0 || mD(l) - dm*s < 0
            continue
        else
            if s > 0 %mk recibe dinero y paga iva
                mD(j) = mD(j) + dm*(1-lambda);
                mD(l) = mD(l) - dm;
                iva = iva + dm*lambda;
            else %ml recibe dinero y paga iva
                mD(j) = mD(j) - dm;
                mD(l) = mD(l) + dm*(1-lambda);
                iva = iva + dm*lambda;
            end
        end
    end
    
    nkD(t,:) = histcounts(mD, Class);
    Sdelta(t) = N*log(N) - sum(nkD(t,:).*log(nkD(t,:)), "omitnan");
    OcaD(t) = sum(a*myCenters.*nkD(t,:));
    OcbD(t) = sum(nkD(t,:).*(1-exp(-a*myCenters)));
    
    VEm = 0;
    for u = 1:C
        VEm = VEm + myCenters(u)*nkD(t,u)/N;
    end
    B(t) = VEm;
    VEm = 0;
end
% 
% hC = histogram(incomeC,CBin1, "Normalization", "cdf");
% cdfC = hC.Values;
% cdfC_vec = 0.5*(hC.BinEdges(1:n-1) + hC.BinEdges(2:n)); 

%% Distribucion Uniforme

iva = 0;
mU = linspace(0, mmax, N); %Distribucion Uniforme
mU = mU/(sum(mU)/M);
Sunif = zeros(1,T);

for t=1:T
    %{
    histogram(mU, Class)
    grid on
    ylim([0 150])
    set(gca, "FontSize", 13, "FontName", "Times")
    xlabel("Dinero (\$)", "Interpreter", "latex")
    ylabel("N\'umero de personas ($n_k$)", "Interpreter", "latex")
    text = sprintf("Distribuci\'on del dinero, t = %d", t);
    title(text, "Interpreter", "latex")
    %title("Distribuci\'on del dinero", "Interpreter", "latex")
    pause(1/120)
    %}
    if mod(t,tau) == 0
        mU = mU + iva/N;
        iva = 0;
    end
    for j=1:N
        l = randi([1 N],1,1);
        s = sign(2*rand()-1);
        if mU(j) + dm*s < 0 || mU(l) - dm*s < 0
            continue
        else
            if s > 0 %mk recibe dinero y paga iva
                mU(j) = mU(j) + dm*(1-lambda);
                mU(l) = mU(l) - dm;
            else %ml recibe dinero y paga iva
                mU(j) = mU(j) - dm;
                mU(l) = mU(l) + dm*(1-lambda);
            end
            iva = iva + dm*lambda;
        end
    end
    nkU = histcounts(mU, Class);
    Sunif(t) = N*log(N) - sum(nkU.*log(nkU), "omitnan");
    OcaU(t) = sum(a*myCenters.*nkU);
    OcbU(t) = sum(nkU.*(1-exp(-a*myCenters)));
end

%{
% ---------- Gráficas Entropía -----------
figure()
plot(Sdelta, "LineWidth", 1.3)
hold on
plot(Sunif, "Linewidth", 1.3)
grid on
axis square
set(gca, "FontSize", 13, "FontName", "Times")   
xlabel("Tiempo $t$", "Interpreter", "latex")
ylabel("Entrop\'ia $S(t)$", "Interpreter", "latex")
%text = sprintf("Reconstruccion a = %d", b);
title("Entrop\'ia en funci\'on del tiempo", "Interpreter", "latex")
%legend("Delta", "Uniforme", "Interpreter", "latex")
%}

%% Graficas Bienestar
%{
figure()
plot(OcbD,"Linewidth", 1.5)
hold on
plot(OcbU,"Linewidth", 1.5)
grid on
axis square
set(gca, "FontSize", 13, "FontName", "Times")
xlabel("Tiempo $t$", "Interpreter", "latex")
ylabel("Bienestar Colectivo $O_2(n_k)$", "Interpreter", "latex")
title("Funcion Objetivo", "Interpreter", "latex")
legend("Delta, $b = 1/100$", "Uniforme, $b = 1/100$", "Interpreter", "latex")
%}

%% Graficas Entropia
%{
figure()
plot(Sunif, "LineWidth", 1.3)
ylim([0, 3200])
hold on
grid on
axis square
set(gca, "FontSize", 13, "FontName", "Times")
xlabel("Tiempo $t$", "Interpreter", "latex")
ylabel("Entrop\'ia $S(t)$", "Interpreter", "latex")
text = sprintf("Delta m = %d", dm);
title("Entrop\'ia en funci\'on del tiempo", "Interpreter", "latex")
legend(text, "Interpreter", "latex")
%}

% Delta m = M*4e-5
% Delta m = M*1e-5
% Delta m = M*8e-5

%% Funcion objetivo con distribucion Boltzmann-Gibbs
%
% ----------- b constante ------------
num = 200;
a = 1;
Temp = 100;
b = 1/Temp;
a_vec = linspace(0, a,num);
O1_bCte = a_vec./b;
O2_bCte = (a_vec./(a_vec+b));
%{
figure()
plot(a_vec, O1_bCte, "LineWidth", 1.8)
hold on
plot(a_vec, O2_bCte, "LineWidth", 1.8)
grid on
set(gca, "FontSize", 13, "FontName", "Times")   
title("Bienestar Colectivo", "Interpreter", "latex")
xlabel("a", "Interpreter", "latex")
ylabel("$O(M)$", "Interpreter", "latex")
legend("Caso A", "Caso B", "Interpreter", "latex")
%}
% ----------- a constante ---------
%a = 1/100;
% 1/200, 1/100;
b_vec = linspace(0, b, num);
O2_aCte = (a./(a+b_vec));
%{
figure()
plot(a_vec, O2_aCte, "LineWidth", 1.8)
grid on
set(gca, "FontSize", 13, "FontName", "Times")   
title("Bienestar Colectivo", "Interpreter", "latex")
xlabel("b", "Interpreter", "latex")
ylabel("$O(M)$", "Interpreter", "latex")
%}

% ----------- a/b y b/a --------------
 a = 1/300;
 b = 1/100;
aeb_vec = linspace(0, b/a, num);
bea_vec = linspace(0, b/a, num);
aEb = aeb_vec./(aeb_vec + 1);
bEa = 1./(1 + bea_vec); 
aEb_O1 = aeb_vec;
 
%{
figure()
plot(aeb_vec, aEb, "LineWidth", 1.8)
hold on
plot(bea_vec, bEa, "LineWidth", 1.8)
plot(aeb_vec, aEb_O1, "LineWidth", 1.8)
% ylim([0, 1])
% xlim([0 b/a])
grid on
set(gca, "FontSize", 13, "FontName", "Times")   
title("Bienestar Colectivo", "Interpreter", "latex")
xlabel("a/b o b/a", "Interpreter", "latex")
ylabel("Bienestar colectivo $O(M)$", "Interpreter", "latex")
title("Boltzmann-Gibbs b/a y a/b")
legend("$O_2(a/b)$", "$O_2(b/a)$", "$O_1(a/b)$", "Interpreter", "latex")
%}