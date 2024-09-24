Q = 0:1:100;
P_a_Q = zeros(5, length(Q));
t = 20:20:100;
mu = sqrt(t);

for t = 1:5
    P_a_Q(t,:) = exp(-mu(t)).*mu(t).^Q/factorial(Q);
end

length(t)
length(P_a_Q)

plot(Q,P_a_Q)
