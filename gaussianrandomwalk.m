close all
x_0 = 0;
x = [x_0];
dt = 0.01;
t = 0:dt:10;
D = 1;
for i =2:length(t)
    a = randn;
    x(i) =x_0+a*sqrt(2*D*dt);
    x_0 = x_0 + a*sqrt(dt);
end
x;
length(x);
length(t);
fig1=figure;
clf;
plot(x,t)
grid on
%axis([-100 100 0 100])
title('Gaussian Random Walk', 'Fontsize', 17, 'Interpreter', 'latex')
xlabel('$x$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$t$', 'Fontsize', 20, 'Interpreter', 'latex')
 
%hold on


