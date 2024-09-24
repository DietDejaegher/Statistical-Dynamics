close all
x = -10:0.01:10;
y = -5:2:-1;
t = 3:5;

%dist = normpdf(x,y,sqrt(t));
clf;
hold on
plot(x, normpdf(x,0,sqrt(1)))
plot(x, normpdf(x,0,sqrt(3)))
for i= 1:length(y)
    plot(x, normpdf(x,y(i),t(i)))
end
ylim([0 0.5])
grid on
xlabel('$x$','Interpreter', 'latex', 'Fontsize', 20)
ylabel('$P(x,y,t)$','Interpreter', 'latex', 'Fontsize', 20)
title('Probability density for diffusive particles', 'Interpreter', 'latex', 'Fontsize', 17)
legend('$y=0\qquad\sigma^2=1$','$y=0\qquad\sigma^2=3$','$y=-5\qquad\sigma^2=3$','$y=-3\qquad\sigma^2=4$','$y=-1\qquad\sigma^2=5$', 'Interpreter', 'latex', 'Location', 'NorthEast', 'Fontsize', 12)