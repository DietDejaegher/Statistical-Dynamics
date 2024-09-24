N = 100;
n = 1:N;
L = 100;
dt = 0.1;
t = 0:dt:1000; 
length(t);%=101
%Matrix: each column corresponds to 1 particle
%each row corresponds to a position of that particle
%The first row are the initial positions of all the particles
%For exemple: x(i,j), i-th row and j-th column, corresponds to
%position x(j) of the i-th particle
x = zeros(length(t),N);
size(x);%101 rows 100 columns
x(1,1) = -N/L; %initial position of 1st particle
size(x(1,1));%scalar 1x1
x(1,1); %=1
for i = 2:N
    x(1,i) = x(1,i-1)-N/L;
end
x(1,:);%intial positions of all the particles
size(x(1,:)); %1 row 100 columns

%Diffusion
for k = 1:N
    for j = 2: length(t)
        a = randn;
        x(j,k) = x(j-1,k) + a*sqrt(dt);
    end
end
%plot(t, x)
%grid on
%xlabel('$t$', 'Fontsize', 20, 'Interpreter', 'latex')
%ylabel('$x$', 'Fontsize', 20, 'Interpreter', 'latex')

%Count the flux, in other words, count the number of particles that ended
%up on the positive side of the origin
x(101,:); %final positions of all the particles
Q = 0;
for h = 1:N
    if x(length(t),h)< 0 
        Q = Q;
    else
        Q = Q+1;
    end
    Q;
end
flux = Q;

%Now we calculate the flux Q a numerous amount of times and for each
%recieved value we give a certain interval a 'point'
%For exemple: The received flux Q=16 we give the interval [16,17[ a point
amount = 10000;
Fluxlist=zeros(1,amount);
for g = 1:length(Fluxlist)
    for k = 1:N
        for j = 2: length(t)
            a = randn;
            x(j,k) = x(j-1,k) + a*sqrt(dt);
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
Fluxlist;

point_per_int = zeros(1,N); 

for i = 1:N
    point_per_int(i) = sum(Fluxlist(:)==i);
end

point_per_int; %the amount of times the flux took the value of each integer between 0 and N

bar(n,point_per_int)
xlabel('$Q$', 'Fontsize', 20, 'Interpreter', 'latex')
title('\textbf{Distribution of diffusive flux}','Fontsize', 17, 'Interpreter', 'latex')


