####################################################
# Gillespie simulations of a simple Poisson Process
####################################################

using LinearAlgebra
using Random
using Pkg
Pkg.add("Plots")
using Plots

N = 10000;

x = zeros(N);
t = zeros(N);
λ = 10;
τ = 1;
x[1] = round(λ * τ);
for i = 2 : N
    reaction_rate = [λ; x[i-1]/τ]
    u_t = rand();
    waiting_time = log(1/u_t) / sum(reaction_rate);
    t[i] = t[i-1] + waiting_time;

    cum_reac = cumsum(reaction_rate) / sum(reaction_rate);
    u_r = rand();
    if u_r <= cum_reac[1]
        x[i] = x[i-1] + 1;
        continue
    elseif u_r <= cum_reac[2]
        x[i] = x[i-1] - 1;
    end
end

plot(t, x, linetype=:steppost)
