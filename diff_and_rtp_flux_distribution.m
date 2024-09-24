close all
tic;
N = 100;
n = 1:N;
L = N;
D = 1;
dt = 0.1;
alpha = 0.5;
prob_switch = alpha*dt; % 0< prob_switch <1
v_0 = 1;
T = 1000;
t = 0:dt:T; 
amount = 1000;

Diffusion_Fluxlist=zeros(1,amount);
RPT_Fluxlist=zeros(1,amount);

for g = 1:amount
    x = zeros(1,N);
    x_diff = zeros(length(t),N);
    x_rtp = zeros(length(t),N);
    for i = 1:N
    position = unifrnd(-L,0);
    x(1,i) = position;
    end
    x_diff(1,:) = x;
    x_rtp(1,:) = x;
    for k = 1:N
        b = rand;
        if b<0.5
            sigma = +1;
        else
            sigma = -1;
        end
        for j = 2: length(t)
            a = randn;
            x_diff(j,k) = x_diff(j-1,k) + a*sqrt(2*D*dt);
            x_rtp(j,k) = x_rtp(j-1,k) + sigma*v_0*dt;
            c = rand;
            if c<prob_switch
                sigma = -sigma;
            else
                sigma = sigma;
            end
        end
    end
    Q_diff = 0;
    for h = 1:N
        if x_diff(length(t),h)< 0 
            Q_diff = Q_diff;
        else
            Q_diff = Q_diff+1;
        end
    end
    Diffusion_Fluxlist(g) = Q_diff;
    
    Q_rtp = 0;
    for h = 1:N
        if x_rtp(length(t),h)<0
            Q_rtp = Q_rtp;
        else
            Q_rtp = Q_rtp+1;
        end
    end
    RTP_Fluxlist(g) = Q_rtp;
end

point_per_int_diff = zeros(1,N);
point_per_int_rtp = zeros(1,N);

for i = 1:N
    point_per_int_diff(i) = sum(Diffusion_Fluxlist(:)==i);
    point_per_int_rtp(i) = sum(RTP_Fluxlist(:)==i);

end

diff_flux_mean = mean(Diffusion_Fluxlist)
rtp_flux_mean = mean(RTP_Fluxlist)

diff_analytical_mean = (N/L)*sqrt(D*T/pi);
rtp_analytical_mean = (N/L)*(v_0/2)*T*exp(-alpha*T)*(besseli(0,alpha*T)+besseli(1,alpha*T));

diff_poiss = zeros(1,length(n));
diff_poiss_analytical=zeros(size(diff_poiss));
rtp_poiss = zeros(1,length(n));
rtp_poiss_analytical = zeros(size(rtp_poiss));

for i = 1:length(n)
    diff_poiss(i) = diff_flux_mean.^i/factorial(i)*exp(-diff_flux_mean);
    diff_poiss_analytical(i) = diff_analytical_mean.^i/factorial(i)*exp(-diff_analytical_mean);
    rtp_poiss(i) = rtp_flux_mean.^i/factorial(i)*exp(-rtp_flux_mean);
    rtp_poiss_analytical(i) = rtp_analytical_mean.^i/factorial(i)*exp(-rtp_analytical_mean);
end

poiss = zeros(2,length(n));

fig1 = figure;
clf;
hold on
for i = 1:length(n)
    plot(n(i), point_per_int_diff(i)/amount, '*b', n(i),point_per_int_rtp(i)/amount, 'xr')
end
plot(n, [diff_poiss_analytical;rtp_poiss_analytical])
%plot(n, rtp_poiss_analytical, 'r')
xlim([0 N])
grid on
title('Probability distribution of the flux', 'Interpreter', 'latex', 'Fontsize', 17)
ylabel('$P_a(Q)$', 'Interpreter', 'latex', 'Fontsize', 20)
xlabel('$Q$','Interpreter', 'latex', 'Fontsize', 20)
legend('Diffusion', 'RTP','Location','NorthEast','Interpreter', 'latex', 'Fontsize', 12) 

fig2 = figure;
clf;
for i = 1:length(n)
    plot(n(i), point_per_int_diff(i)/amount, '*b', n(i),point_per_int_rtp(i)/amount, 'xr')
    hold on
end
plot(n, [diff_poiss_analytical;rtp_poiss_analytical])
xlim([0 35])
grid on
title('Probability distribution of the flux', 'Interpreter', 'latex', 'Fontsize', 17)
ylabel('$P_a(Q)$', 'Interpreter', 'latex', 'Fontsize', 20)
xlabel('$Q$','Interpreter', 'latex', 'Fontsize', 20)
set(gca, 'YScale', 'log')
%set(gca, 'xtick',0:10:length(n))
legend('Diffusion', 'RTP','Location','NorthEast','Interpreter', 'latex', 'Fontsize', 12) 

fig3=figure;
clf;
for i = 1:length(n)
    plot(n(i), point_per_int_diff(i)/amount, '*b')
    hold on
end
plot(n, diff_poiss_analytical, 'b')
xlim([0 N])
grid on
title('Probability distribution of diffusive flux', 'Interpreter', 'latex', 'Fontsize', 17)
ylabel('$P_a(Q)$', 'Interpreter', 'latex', 'Fontsize', 20)
xlabel('$Q$','Interpreter', 'latex', 'Fontsize', 20)
%legend('Diffusion','Location','NorthEast','Interpreter', 'latex', 'Fontsize', 12) 


fig4=figure;
clf;
for i = 1:length(n)
    plot(n(i),point_per_int_rtp(i)/amount, 'xr')
    hold on
end
plot(n,rtp_poiss_analytical,'r')
xlim([0 N])
grid on
title('Probability distribution of RTP flux', 'Interpreter', 'latex', 'Fontsize', 17)
ylabel('$P_a(Q)$', 'Interpreter', 'latex', 'Fontsize', 20)
xlabel('$Q$','Interpreter', 'latex', 'Fontsize', 20)
%legend('RTP','Location','NorthEast','Interpreter', 'latex', 'Fontsize', 12) 

fig5=figure;
clf;
for i = 1:length(n)
    plot(n(i), point_per_int_diff(i)/amount, '*b')
    hold on
end
plot(n, diff_poiss_analytical, 'b')
xlim([0 35])
grid on
title({'Probability distribution of diffusive flux';'with logarithmic $y$-axis'}, 'Interpreter', 'latex', 'Fontsize', 17)
ylabel('$P_a(Q)$', 'Interpreter', 'latex', 'Fontsize', 20)
xlabel('$Q$','Interpreter', 'latex', 'Fontsize', 20)
set(gca,'YScale', 'log')
%legend('Diffusion','Location','NorthEast','Interpreter', 'latex', 'Fontsize', 12) 

fig6=figure;
clf;
for i = 1:length(n)
    plot(n(i),point_per_int_rtp(i)/amount, '-xr')
    hold on
end
plot(n,rtp_poiss_analytical,'r')
xlim([0 35])
grid on
title({'Probability distribution of RTP flux';'with logarithmic $y$-axis'}, 'Interpreter', 'latex', 'Fontsize', 17)
ylabel('$P_a(Q)$', 'Interpreter', 'latex', 'Fontsize', 20)
xlabel('$Q$','Interpreter', 'latex', 'Fontsize', 20)
set(gca,'YScale', 'log')
%legend('RTP','Location','NorthEast','Interpreter', 'latex', 'Fontsize', 12) 
toc;