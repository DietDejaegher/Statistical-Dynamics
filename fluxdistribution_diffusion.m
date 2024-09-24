close all
tic;
N = 100;
N1 = 10;
n = 1:N;
L = N;
L1 =N1;
n1=1:N1;
D = 1;
amount = 1000;
dt = 0.1;
T = 1000;
t = 0:dt:T; 
length(t);%=101
%Matrix: each column corresponds to 1 particle
%each row corresponds to a position of that particle
%The first row are the initial positions of all the particles
%For exemple: x(i,j), i-th row and j-th column, corresponds to
%position x(j) of the i-th particle
x = zeros(length(t),N);
size(x);%101 rows 100 columns
size(x(1,1));%scalar 1x1
for i = 1:N
    position = unifrnd(-L,0);
    x(1,i) = position;
end
x(1,:);%intial positions of all the particles
size(x(1,:)); %1 row 100 columns

%Diffusion
for k = 1:N
    for j = 2: length(t)
        a = randn;
        x(j,k) = x(j-1,k) + a*sqrt(2*D*dt);
    end
end
%plot(t, x)
%grid on
%xlabel('$t$', 'Fontsize', 20, 'Interpreter', 'latex')
%ylabel('$x$', 'Fontsize', 20, 'Interpreter', 'latex')

%Count the flux, in other words, count the number of particles that ended
%up on the positive side of the origin
x(length(t),:); %final positions of all the particles
Q = 0;
for h = 1:N
    if x(length(t),h)< 0 
        Q = Q;
    else
        Q = Q+1;
    end
end
flux = Q;

%Now we calculate the flux Q a numerous amount of times and for each
%recieved value we give a certain interval a 'point'
%For exemple: The received flux Q=16 we give the interval [16,17[ a point



Fluxlist=zeros(1,amount);
for g = 1:amount
    x = zeros(length(t),N);
    for i = 1:N
    position = unifrnd(-L,0);
    x(1,i) = position;
    end
    for k = 1:N
        for j = 2: length(t)
            a = randn;
            x(j,k) = x(j-1,k) + a*sqrt(2*D*dt);
        end
    end
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
Fluxlist; %all the received flux values
mean_flux = mean(Fluxlist)
var_flux = var(Fluxlist)

expected_Q_mean = (N/L)*sqrt(D*T/(pi))
%expected_Q_mean2 = N/L * sqrt(mean(dx2)/(2*pi))

point_per_int = zeros(1,N); 

for i = 1:N
    point_per_int(i) = sum(Fluxlist(:)==i);
end

diff_analytical_mean = (N/L)*sqrt(D*T/pi);

poiss = zeros(1,length(n));
diff_poiss_analytical=zeros(size(poiss));

for i = 1:length(n)
    poiss(i) = mean_flux.^i/factorial(i)*exp(-mean_flux);
    diff_poiss_analytical(i) = diff_analytical_mean.^i/factorial(i)*exp(-diff_analytical_mean);
end
point_per_int; %the amount of times the flux took the value of each integer between 0 and N


fig1 = figure;
clf;
bar(n,point_per_int/amount)
hold on
plot(n, poiss)
xlim([0 N])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title('Distribution of diffusive flux','Fontsize', 17, 'Interpreter', 'latex')


fig2 = figure;
clf;
bar(n,point_per_int/amount)
hold on
plot(n, poiss)
set(gca, 'YScale', 'log')
%ylim([0 1])
xlim([0 N])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Distribution of diffusive flux';'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')

fig3 = figure;
hold on
grid on
plot(n, poiss)
set(gca, 'YScale', 'log')
xlim([0 N])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Distribution of diffusive flux';'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')

fig4= figure;
clf;
hold on
grid on
for i = 1:length(n)
    plot(n(i), point_per_int(i)/amount, '*b')
end
plot(n, diff_poiss_analytical, 'b')
xlim([0 N])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title('Distribution of diffusive flux','Fontsize', 17, 'Interpreter', 'latex')



fig5 = figure;
clf;
hold on
grid on
for i = 1:length(n)
    plot(n(i), point_per_int(i)/amount, '*b')
end
plot(n, diff_poiss_analytical, 'b')
xlim([0 40])
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('$P_a(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Distribution of diffusive flux';'with logarithmic $y$-axis'},'Fontsize', 17, 'Interpreter', 'latex')
set(gca, 'YScale', 'log')

toc;