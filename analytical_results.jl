import Pkg
Pkg.activate("packages")
Pkg.instantiate()
using Roots, SpecialFunctions

critical_number_mutations(w0, s) = -log(w0)/log(1+s)

reproduction_rate(m, w0, s) = w0 * (1+s)^m

population_size_meltdown(t, w0, s, mu, K) = K*(1+s)^(mu/2 * t * (t-1))

mean_num_mutations(t, s, mu) = -mu * (1+s)/s * (1 - (1+s)^t)

mutation_selection_balance(s, mu) = -mu*(1+s)/s

maximal_s(mu, K) = mu/(log(K) - mu)

# Muller's ratchet

function start_ratchet(w0, s, mu, K)
    if K*exp(mu*(1+s)/s) < 1
        return log(1 + s*log(K)/(mu*(1+s)))/log(1+s)
    else
        return Inf # the ratchet never starts
    end
end
function start_ratchet(w0, s, mu, K, discrete)
    if K*exp(mu*(1+s)/s) < 1
        return round(Int, log(1 + s*log(K)/(mu*(1+s)))/log(1+s), RoundUp)
    else
        return Inf # the ratchet never starts
    end
end

speed_ratchet_Lynch(s, mu, K) = mu*(1+s)/(1-K*s)

speed_ratchet_Rouzine(s, mu, K) = 1 + s*log(K)/mu

function speed_ratchet_Pamilo(s, mu, K)
    if mu > -s*log(K)
        return mu + s*log(K)
    else
        return 0
    end
end

function speed_ratchet_Gessler(s, mu, K)
    k = 0
    while log(K) + mu*(1+s)/s + k*log(-mu*(1+s)/s) - log(factorial(big(k))) < 0
        k += 1
    end
    b = k
    A = log(K) + (b-k)*log(-mu*(1+s)/s - k) - log(factorial(big(b - k))) - (-mu*(1+s)/s - k)
    while A < 0
        b += 1
        a = 0
        for i = 0:b-k-1
            a += (-mu*(1+s)/s-k)^i / factorial(big(i))
        end
        A = log(K) + (b-k)*log(-mu*(1+s)/s - k) - log(factorial(big(b - k))) -  log(exp(-mu*(1+s)/s - k) - a)
    end
    return -s*(2*k - b)
end

function speed_ratchet_Gessler(s, mu, K, continuous)
    global mu_dummie = mu
    global s_dummie = s
    global K_dummie = K
    f_k(x) = log(K_dummie) + mu_dummie*(1+s_dummie)/s_dummie + x*log(-mu_dummie*(1+s_dummie)/s_dummie) - loggamma(x+1)
    function f_b(x)
        S = 0
        for i = 0:Int(round(x-k-1, RoundDown))
            S += (-mu_dummie*(1+s_dummie)/s_dummie - k)^i / factorial(big(i))
        end
        return log(K_dummie) + (x-k)*log(-mu_dummie*(1+s_dummie)/s_dummie - k) - loggamma(x-k+1) - log(exp(-mu_dummie*(1+s_dummie)/s_dummie - k) - S - ((x-k-1) - Int(round(x-k-1, RoundDown))) * (-mu_dummie*(1+s_dummie)/s_dummie - k)^Int(round(x-k, RoundDown)) / factorial(big(Int(round(x-k, RoundDown)))))
    end
    try
        global k = find_zero(f_k, (0, -mu_dummie*(1+s_dummie)/s_dummie))
    catch e
        global k = 0
    end
    try
        global b = find_zero(f_b, (k, -mu_dummie*(1+s_dummie)/s_dummie))
    catch e
        global b = 0
    end
    return -s_dummie*(2*k - b)
end

# (Extinction) time

ratchet_time(w0, s, mu, K) = (-log(w0)/log(1+s) - log(K))/speed_ratchet_Gessler(s, mu, K, 1)

meltdown_time(s, mu, K) = (-2*log(K)/(mu * log(1+s)) + 1/4)^0.5 - 1/2

extinction_time(w0, s, mu, K) = start_ratchet(w0, s, mu, K) + (-log(w0)/log(1+s) - log(K))/speed_ratchet_Gessler(s, mu, K, 1) + (-2*log(K)/(mu * log(1+s)) + 1/4)^0.5 - 1/2

function optimal_s(w0, mu, K, low_s, high_s)
    global w0_dummie = w0
    global mu_dummie = mu
    global K_dummie = K
    extinction_time_s(s) = extinction_time(w0_dummie, s, mu_dummie, K_dummie)
    res = optimize(extinction_time_s, low_s, high_s)
    return Optim.minimizer(res)
end
