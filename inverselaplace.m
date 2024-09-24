x_0 = 0;
x = [x_0];
t = 0:0.01:10;
for i =2:length(t)
    a = randn;
    x(i) =x_0+a;
    x_0 = x_0 + a;
end
x;
length(x)
length(t)

plot(x,t)
xlabel('$x$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$t$', 'Fontsize', 20, 'Interpreter', 'latex')



