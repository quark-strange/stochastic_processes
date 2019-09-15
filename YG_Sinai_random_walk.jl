####################################################
# Y.G. Sinai's random walk in random field
####################################################

using LinearAlgebra
using Random
using Plots

# Generate a Landscape
function landscape_generator()
    x_max = 5 * 10^4;
    δ = 2;
    U = rand([δ, -δ],x_max);
    return cumsum(U);
end

# Random Walker
function random_walker(T, U)
    x = ones(T);
    x[1] = length(U)/2;
    x = convert(Array{Int64}, x);
    for i = 2 : T
        δx = rand([1, -1]);
        ΔU = U[x[i-1] + δx] - U[x[i-1]];
        if ΔU <= 0 || rand() <= exp(-ΔU)
            x[i] = x[i-1] + δx;
        else
            x[i] = x[i-1];
        end
    end
    return (x - x[1]*ones(T)).^2
end

# Generate a single MSD
function single_MSD(T, ens)
    x_ms = zeros(T);
    for i = 1 : ens
        x_ms;
        x_ms += random_walker(T, landscape_generator()) ./ ens;
    end
    return x_ms
end

let
    instance_Potential = 300;
    T = 5 * 10^4;
    t = [n for n in 1:T];
    ens = 100;
    MSD_s = zeros(T); for i = 1 : instance_Potential
        MSD_s += single_MSD(T, ens) / instance_Potential;
    end
    global MSD_s;
    global t;
    plot((log.(t)).^4, MSD_s, legend=:none)
end
