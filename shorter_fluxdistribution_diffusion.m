close all
tic;
N = 100;
N1 =2000;
n1 = 1:N1;
L1=N1;
n = 1:N;
L = N;
D = 1;
amount = 100;
dt = 0.1;
T = 1000;
t = 0:dt:T; 

x = zeros(length(t),N);
x1= zeros(length(t),N1);
Fluxlist=zeros(1,amount);
Fluxlist1=zeros(1,amount);


for g = 1:amount
    x = zeros(1,N);
    x1= zeros(length(t),N1);
    for i = 1:N
    position = unifrnd(-L,0);
    x(1,i) = position;
    end
    
    for i = 1:N1
    position1 = unifrnd(-L1,0);
    x1(1,i) = position1;
    end
    
    for k = 1:N
        for j = 2: length(t)
            a = randn;
            x(1,k) = x(1,k) + a*sqrt(2*D*dt);
        end
    end
    
    for k = 1:N1
        for j = 2: length(t)
            a1 = randn;
            x1(1,k) = x1(1,k) + a1*sqrt(2*D*dt);
        end
    end
    
    Q = 0;
    Q1=0;
    for h= 1:N
        if x(1,h) <0
            Q=Q;
        else
            Q=Q+1;
        end
    end
    
    for h= 1:N1
        if x1(1,h) <0
            Q1=Q1;
        else
            Q1=Q1+1;
        end
    end
    Fluxlist(g)=Q;
    Fluxlist1(g)=Q1;
end

mean_flux = mean(Fluxlist)
mean_flux1 = mean(Fluxlist1)

var_flux=var(Fluxlist)
var_flux1=var(Fluxlist1)


expected_Q_mean = (N/L)*sqrt(D*T/pi)

point_per_int = zeros(1,N);
point_per_int1 = zeros(1,N1);


for i = 1:N
    point_per_int(i) = sum(Fluxlist(:)==i);
end

for i = 1:N1
    point_per_int1(i) = sum(Fluxlist1(:)==i);
end


poiss = zeros(1,length(n));
poiss1 = zeros(1,length(n1));
poiss_analytical = zeros(1, length(n));
for i = 1:length(n)
    poiss(i) = mean_flux.^i/factorial(i)*exp(-mean_flux);
    
end
point_per_int; %the amount of times the flux took the value of each integer between 0 and N

for i = 1:length(n1)
    poiss1(i) = mean_flux1.^i/factorial(i)*exp(-mean_flux1);
    poiss_analytical(i) = exp(-expected_Q_mean)*expected_Q_mean^i/factorial(i);
end

fig1 = figure;
clf;
bar(n,point_per_int/amount)
hold on
plot(n, poiss)
xlim([0 100])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title('Distribution of diffusive flux','Fontsize', 17, 'Interpreter', 'latex')


fig2 = figure;
clf;
bar(n,point_per_int/amount)
hold on
plot(n, poiss)
set(gca, 'YScale', 'log')
%ylim([0 1])
xlim([0 100])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Distribution of diffusive flux';'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')

fig3 = figure;
hold on
grid on
plot(n, poiss, '*b')
plot(n1, poiss1, '+r')
plot(n1, poiss_analytical)
set(gca, 'YScale', 'log')
xlim([0 500])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Distribution of diffusive flux';'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')

toc;



    
    
    