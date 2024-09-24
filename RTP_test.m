x_0 = 0;
a = 0.5;%average rate of flipping
v_0 = 1;%constant velocity
x = [x_0];
dt = 0.01;
t = 0:dt:10;
b = rand;
if b<0.5
    sigma = +1;
else
    sigma = -1;
end

for i = 2:length(t)
    x(i) = x_0 + sigma*v_0*dt;
    x_0 = x_0 + sigma*v_0*dt;
    c = randn;
    if c<a*dt*(1-a*dt)
        sigma = -sigma;
    else
        sigma = sigma;
    end
end

plot(x,t)
grid on
title('Run-And-Tumble Walk', 'Fontsize', 17, 'Interpreter', 'latex')
xlabel('$x$',  'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$t$', 'Fontsize', 20, 'Interpreter', 'latex')