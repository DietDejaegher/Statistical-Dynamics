close all
a = 0.5;
v_0 = sqrt(2*a);
t = 0:0.01:200;
D_eff = v_0^2/(2*a);
I = zeros(3,length(t));
for nu = 0:1
    I(nu+1,:) = besseli(nu,a*t);
end
I(3,:) = I(1,:)+I(2,:);

fig1 = figure ;
clf;
plot(t,I)
legend('$I_0$','$I_1$', '$I_0+I_1$','Location','NorthWest','Interpreter', 'latex', 'Fontsize', 12)
grid on
title('Modified Bessel Functions of the First Kind','Fontsize', 17, 'Interpreter', 'latex')
xlabel('$t$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$I_n(\alpha t)$', 'Fontsize', 20, 'Interpreter', 'latex')
%hold on

size(I(3,:));
size(exp(t));
mu = (v_0/2)*t.*exp(-a*t).*I(3,:);

fig2 = figure ;
plot(t, mu)
grid on
title('$\mu(t)$ for RTP', 'Fontsize', 17, 'Interpreter', 'latex')
xlabel('$t$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$\mu(t)$', 'Fontsize', 20, 'Interpreter', 'latex')

muoft = zeros(2,length(t));
muoft(1,:)=mu;
muoft(2,:)= sqrt(D_eff*t/pi);
muoft_rtp= muoft(1,:);

fig3 = figure;
clf;
hold on
plot(t, muoft)
%plot(t, v_0/2*t.*(1-a*t/2), '--r')
grid on
xlim([100 200])
xlabel('$t$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$\mu(t)$','Fontsize', 20, 'Interpreter', 'latex')
title('$\mu(t)$ for diffusion and RTP','Fontsize', 17, 'Interpreter', 'latex')
legend('RTP', 'Diffusion', 'Location', 'NorthWest','Fontsize', 12, 'Interpreter', 'latex')

fig4 = figure;
plot(t, muoft(2,:))
grid on
xlabel('$t$','Fontsize', 20, 'Interpreter', 'latex')
ylabel('$\mu(t)$','Fontsize', 20, 'Interpreter', 'latex')
title('$\mu(t)$ for diffusion','Fontsize', 17, 'Interpreter', 'latex')
