close all
syms k x
q=[0.2 0.5 1 2 ] ;
D = 1;
a = 0.5;
v_0 = sqrt(2*a);
rho = 1;
flux_val=0:1:100;
t0 = 0:0.1:1000;
t1 = 0:1:1000;
time = [10 100 500 1000];
Q_small = [0 5 10 20];
Q_large = [100 200 500 1000];
signs = ['*' 'x' '+' 'o'];
colors = ['b' 'r' 'g' 'm'];
name = {'t=10', 't=100', 't=500', 't=1000'};


mu_time1 = sqrt(D*(time)/pi);
mu_time2 = rho*v_0/2*(time).*exp(-a*time).*(besseli(0,a*time)+besseli(1,a*time));

mu_t_diff = sqrt(D*(t0)/pi);
mu_t_rtp = rho*v_0/2*(t0).*exp(-a*t0).*(besseli(0,a*t0)+besseli(1,a*t0));
mu_t_rtp_discrete = rho*v_0/2*(t1).*exp(-a*t1).*(besseli(0,a*t1)+besseli(1,a*t1));


P_Q_diff = zeros(length(time),length(flux_val));%Probability distribution in function of the flux Q
P_Q_rtp = zeros(size(P_Q_diff));
P_t_diff = zeros(length(Q_small),length(t0));%Probability distribution in function of time t
P_t_rtp = zeros(size(P_t_diff));
P_large_deviation_diff = zeros(length(q),length(t0));
P_max_rtp = exp(-mu_t_rtp).*mu_t_rtp.^(rho*v_0*t0)./gamma(rho*v_0*t0+1);
P_max_rtp_discrete = exp(-mu_t_rtp_discrete).*mu_t_rtp_discrete.^(rho*v_0*t1)./gamma(rho*v_0*t1+1);

%P_Q_test= exp(-mu_time1)*mu_time1.^flux_val/factorial(flux_val);
%size(P_a);

for j = 1:length(time)
    for i = 1:length(flux_val)
        P_Q_diff(j,i)= exp(-mu_time1(j))*mu_time1(j)^flux_val(i)/factorial(flux_val(i));
        P_Q_rtp(j,i) = exp(-mu_time2(j))*mu_time2(j)^flux_val(i)/factorial(flux_val(i));
    end
end

for j = 1:length(Q_small)
    P_t_diff(j,:) = exp(-mu_t_diff).*mu_t_diff.^Q_small(j)/factorial(Q_small(j));
    P_t_rtp(j,:) = exp(-mu_t_rtp).*mu_t_rtp.^Q_small(j)/factorial(Q_small(j));
end

for j = 1:length(q)
    P_large_deviation_diff(j,:) = exp(-mu_t_diff*(q(j)*log(q(j))-q(j)+1));
    %P_t_rtp(j,:) = exp(-mu_t_rtp).*mu_t_rtp.^Q_small(j)/factorial(Q_small(j));
end
        
fig1 = figure;
clf;
for j = 1:length(time)
    for i = 1:length(flux_val)
        plot(flux_val(i), P_Q_diff(j,i), signs(j), 'color',colors(j))
        hold on
    end
    h(j) = plot(flux_val,P_Q_diff(j,:), colors(j));
    hold on
end
%ylim([0 1])
title({'$P_a(Q)$ for different values of $t$'; 'in case of diffusion'},'Fontsize', 17, 'Interpreter', 'latex')
xlim([0 flux_val(length(flux_val))/2])
xlabel('$Q$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
legend(h,{'$t=10$','$t=100$', '$t=500$', '$t=1000$'}, 'Location', 'NorthEast','Fontsize', 12, 'Interpreter', 'latex') 
grid on



fig2 = figure;
clf;
for j = 1:length(time)
    for i = 1:length(flux_val)
        plot(flux_val(i), P_Q_rtp(j,i), signs(j), 'color',colors(j))
        hold on
    end
    h(j) = plot(flux_val,P_Q_rtp(j,:), colors(j));
    hold on
end
%ylim([0 1])
title({'$P_a(Q)$ for different values of $t$'; 'in case of RTP'},'Fontsize', 17, 'Interpreter', 'latex')
xlim([0 flux_val(length(flux_val))/2])
xlabel('$Q$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
legend(h,{'$t=10$','$t=100$', '$t=500$', '$t=1000$'}, 'Location', 'NorthEast','Fontsize', 12, 'Interpreter', 'latex') 
grid on


fig3 = figure;
clf;
for j = 1:length(Q_small)
    h(j) = plot(t0,P_t_diff(j,:), colors(j));
    hold on
end
%ylim([0 1])
title({'$P_a(t)$ for small values of $Q$'; 'in case of diffusion'},'Fontsize', 17, 'Interpreter', 'latex')
xlim([0 t0(length(t0))])
xlabel('$t$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(t)$', 'Fontsize', 20, 'Interpreter', 'latex')
legend(h,{'$Q=0$','$Q=5$', '$Q=10$', '$Q=20$'}, 'Location', 'NorthEast','Fontsize', 12, 'Interpreter', 'latex') 
grid on

fig4 = figure;
clf;
for j = 1:length(Q_small)
    h(j) = plot(t0,P_t_rtp(j,:), colors(j));
    hold on
end
%ylim([0 1])
title({'$P_a(t)$ for different small of $Q$'; 'in case of RTP'},'Fontsize', 17, 'Interpreter', 'latex')
xlim([0 t0(length(t0))])
xlabel('$t$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(t)$', 'Fontsize', 20, 'Interpreter', 'latex')
legend(h,{'$Q=0$','$Q=5$', '$Q=10$', '$Q=20$'}, 'Location', 'NorthEast','Fontsize', 12, 'Interpreter', 'latex') 
grid on

fig5=figure;
clf;
for j =1:length(Q_large)
    h(j) = plot(t0,P_large_deviation_diff(j,:), colors(j));
    hold on
end
title({'$P_a(t)$ for different large of $Q$'; 'in case of diffusion'},'Fontsize', 17, 'Interpreter', 'latex')
xlim([0 t0(length(t0))])
xlabel('$t$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(t)$', 'Fontsize', 20, 'Interpreter', 'latex')
set(gca, 'YScale', 'log')
legend(h,{'$q=0.2$','$q=0.5$', '$q=1$', '$q=2$'}, 'Location', 'NorthEast','Fontsize', 12, 'Interpreter', 'latex') 
grid on

fig6=figure;
clf;
hold on
h(1)=plot(rho*v_0*t1, P_max_rtp_discrete, '*b');
h(2)=plot(t0, P_max_rtp, 'b');
xlim([0 10])
title('Probability distribution for $Q=Q_{max}$', 'Fontsize', 17, 'Interpreter', 'latex')
xlabel('$t$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(t)$', 'Fontsize', 20, 'Interpreter', 'latex')
grid on