close all
D = 1;
v_0 = 1;
alpha = 0.5;
dt = 0.1;
t = 1000;
q = 0:100; %ratio Q/mu
mu_diff = sqrt(D*t/pi);
mu_rtp = v_0/2*t.*exp(-alpha*t).*(besseli(0,alpha*t)+besseli(1,alpha*t));

P_t = zeros(2, length(q));
P_t(1,:) = exp(-mu_diff.*(q/mu_diff.*log(q/mu_diff)-q/mu_diff+1));
P_t(2,:) = exp(-mu_rtp.*(q/mu_rtp.*log(q/mu_rtp)-q/mu_rtp+1));

LDF = zeros(2,length(q));
LDF(1,:) = q/mu_diff.*log(q/mu_diff)-q/mu_diff+1;
LDF(2,:) = q/mu_rtp.*log(q/mu_rtp)-q/mu_rtp+1;

fig1 = figure;
clf;
grid on
plot(q, P_t)
legend('Diffusion', 'RTP', 'Location', 'NorthEast', 'Interpreter', 'latex','Fontsize', 12)
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_\psi(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title('Large Deviation Probability', 'Fontsize', 17, 'Interpreter', 'latex')


fig2 = figure;
clf;
grid on
plot(q, P_t)
legend('Diffusion', 'RTP', 'Location', 'NorthEast', 'Interpreter', 'latex','Fontsize', 12)
set(gca, 'YScale', 'log')
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_\psi(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Large Deviation Probability'; 'with logarithmic $y$-axis'} , 'Fontsize', 17, 'Interpreter', 'latex')

fig3 = figure;
clf;
grid on
plot(q, LDF)
xlabel('$q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$\psi(q)$','Fontsize', 20, 'Interpreter', 'latex')
title('Large Deviation Function','Fontsize', 17, 'Interpreter', 'latex')
