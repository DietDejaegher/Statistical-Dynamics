close all
N = 100;
n = 1:N;
L = N;
dt = 0.1;
alpha = 0.5;
prob_switch = alpha*dt; % 0< prob_switch <1
v_0 = 1;
T =1000;
t = 0:dt:T;

amount = 1000;
Fluxlist=zeros(1,amount);
for g = 1:length(Fluxlist)
    x = zeros(length(t),N);
    for i = 1:N
    position = unifrnd(-L,0);
    x(1,i) = position;
    end
    for k = 1:N
        b = rand;
        if b<0.5
            sigma = +1;
        else
            sigma = -1;
        end
        for j = 2: length(t)
            x(j,k) = x(j-1,k) + sigma*v_0*dt;
            c = rand;
            if c<prob_switch
                sigma = -sigma;
            else
                sigma = sigma;
            end
        end
    end
    %plot(t,x)
    Q = 0;
    for h = 1:N
        if x(length(t),h)< 0 
            Q = Q;
        else
            Q = Q+1;
        end
    end
    Fluxlist(g) = Q;
    
end
Fluxlist;
mean_flux = mean(Fluxlist)
var_flux = var(Fluxlist)

expected_Q_mean = (N/L)*(v_0/2)*T*exp(-alpha*T)*(besseli(0,alpha*T)+besseli(1,alpha*T))
expected_Q_mean1 = (N/L)*sqrt(v_0^2*T/(2*alpha*pi))
point_per_int = zeros(1,N); 

for i = 1:N
    point_per_int(i) = sum(Fluxlist(:)==i);
end

rtp_analytical_mean = (N/L)*(v_0/2)*T*exp(-alpha*T)*(besseli(0,alpha*T)+besseli(1,alpha*T));

poiss = zeros(1,length(n));
rtp_poiss_analytical = zeros(size(poiss));

for i = 1:length(n)
    poiss(i) = mean_flux.^i/factorial(i)*exp(-mean_flux);
    rtp_poiss_analytical(i) = rtp_analytical_mean.^i/factorial(i)*exp(-rtp_analytical_mean);

end

point_per_int; %the amount of times the flux took the value of each integer between 0 and N


fig1 = figure;
clf;
bar(n,point_per_int/amount)
hold on
plot(n, poiss)
xlim([0 100])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title('Distribution of RTP flux','Fontsize', 17, 'Interpreter', 'latex')

fig2 = figure;
clf;
bar(n,point_per_int/amount)
hold on
plot(n, poiss)
xlim([0 N])
set(gca, 'YScale', 'log')
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Distribution of RTP flux';'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')

fig3 = figure;
grid on
plot(n, poiss)
set(gca, 'YScale', 'log')
xlim([0 N])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Distribution of RTP flux';'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')

fig4 = figure;
clf;
hold on
grid on
for i = 1:length(n)
    plot(n(i), point_per_int(i)/amount, 'xr')
end
plot(n, rtp_poiss_analytical, 'r')
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title('Distribution of RTP flux','Fontsize', 17, 'Interpreter', 'latex')

fig5 = figure;
clf;
hold on
grid on
for i = 1:length(n)
    plot(n(i), point_per_int(i)/amount, 'xr')
end
plot(n, rtp_poiss_analytical, 'r')
xlim([0 40])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
set(gca, 'YScale', 'log')
title({'Distribution of RTP flux'; 'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')



