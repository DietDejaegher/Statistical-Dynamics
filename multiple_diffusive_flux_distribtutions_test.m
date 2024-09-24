close all
tic;
N = 1000;
n = 1:N;
L = N;
dt = 0.1;
T = 1000;
t = 0:dt:T; 
length(t);

D = zeros(1,amount);

amount = 100;
Fluxlist=zeros(1,amount);
clf;
for l = 200:200:T
    x = zeros(l/dt+1,N);
    x(1,1) = -N/L;
    for a = 2:N
        x(:,a) = x(1,a-1) - N/L;
    end
    
    for g = 1:amount
        for k = 1:N
            for j = 2: (l/dt+1)
                a = randn;
                x(j,k) = x(j-1,k) + a*sqrt(dt);
            end
        end
        Q = 0;
        for h = 1:N
            if x(l/dt+1,h)< 0 
                Q = Q;
            else
                Q = Q+1;
            end
        end
        Fluxlist(g) = Q;
    
        %Diffusion coefficient is the average of the 
        %displacement^2 divided by 2*dt
        dx2 = zeros(1,N); %displacements squared of all the particles
        for i = 1:N
            %displacement squared of each particle
            dx2(i)= (x(l/dt+1,i) - x(1,i))^2;
        end
        D(g) = mean(dx2)/(2*l);
    end
    Fluxlist; %all the received flux values
    mean_flux = mean(Fluxlist);
    var_flux = var(Fluxlist);

    D_mean = mean(D);

    expected_Q_mean = (N/L)*sqrt(D_mean*l/(pi));
    %expected_Q_mean2 = N/L * sqrt(mean(dx2)/(2*pi))

    point_per_int = zeros(1,N); 

    for i = 1:N
        point_per_int(i) = sum(Fluxlist(:)==i);
    end

    poiss = zeros(1,length(n));
    for i = 1:length(n)
        poiss(i) = mean_flux.^i/factorial(i)*exp(-mean_flux);
    end
    point_per_int; %the amount of times the flux took the value of each integer between 0 and N

    %bar(n,point_per_int/amount)
    plot(n, poiss)
    hold on
    %ylim([0 1])
    xlim([0 100])
    xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
    ylabel('$P(Q)$', 'Fontsize', 20, 'Interpreter', 'latex')
    set(gca, 'YScale', 'log')
    title('Distribution of diffusive flux at different times','Fontsize', 17, 'Interpreter', 'latex')
end
legend('$t=200$', '$t=400$', '$t=600$', '$t=800$', '$t=1000$', 'Fontsize', 12, 'Interpreter', 'latex', 'Location', 'NorthEast')

toc;