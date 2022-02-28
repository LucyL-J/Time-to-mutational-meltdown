include("simulate_mutational_meltdown.jl")

function ratchet_speeds_s_del(run_s_del::Vector{Float64}; rng=Random.seed!(), runs::Int=1)
  # runs the simulation for different deleterious selection coefficients given by -10^run_s_del
  # creates two output files which can be uniquely identified by the seed, one file with the mean ratchet speed and the other with the parameters
  # can be called with a seed to reproduce a certain series of runs

  parameters = DataFrame(runs=runs, w0=wt_reproduction_rate, mu_del=mu_del, N0=founder_population_size, K=carrying_capacity, genmax=max_generations, SEED=rng.seed)
  ratchet_speeds = DataFrame(runs=1:runs)

  mean_speed = Vector{Float64}(undef, runs)
  for i in eachindex(run_s_del)
    global s_del = -10^run_s_del[i] # be careful, the selection coefficient is given by -10^(input values)
    for j = 1:runs
      mean_speed[j] = 1 / mean(deleteat!(times_until_loss_fittest_class(), 1)) # do not include the loss of the zero-mutations class
    end
    ratchet_speeds[:, "$(run_s_del[i])"] = mean_speed
  end

  CSV.write(string("data/ratchet_speeds_s_del_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), ratchet_speeds)
  CSV.write(string("data/ratchet_speeds_s_del_para_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), parameters)
  return ratchet_speeds, parameters
end

function ratchet_speeds_mu_del(run_mu_del::Vector{Float64}; rng=Random.seed!(), runs::Int=1)
  # runs the simulation for different mutation rates given by 10^run_mu_del
  # creates two output files which can be uniquely identified by the seed, one file with the mean ratchet speed and the other with the parameters
  # can be called with a seed to reproduce a certain series of runs

  parameters = DataFrame(runs=runs, w0=wt_reproduction_rate, s_del=s_del, N0=founder_population_size, K=carrying_capacity, genmax=max_generations, SEED=rng.seed)
  ratchet_speeds = DataFrame(runs=1:runs)

  mean_speed = Vector{Float64}(undef, runs)
  for i in eachindex(run_mu_del)
    global mu_del = 10^run_mu_del[i] # be careful, the mutation rate is given by 10^(input values)
    for j = 1:runs
      mean_speed[j] = 1 / mean( deleteat!(times_until_loss_fittest_class(), 1)) # do not include the loss of the zero-mutations class
    end
    ratchet_speeds[:, "$(run_mu_del[i])"] = mean_speed
  end
  CSV.write(string("data/ratchet_speeds_mu_del_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), ratchet_speeds)
  CSV.write(string("data/ratchet_speeds_mu_del_para_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), parameters)
  return ratchet_speeds, parameters
end

function ratchet_speeds_K(run_K::Vector{Float64}; rng=Random.seed!(), runs::Int=1)
  # runs the simulation for different carrying capacities given by 10^run_K rounded to an Integer
  # creates two output files which can be uniquely identified by the seed, one file with the mean ratchet speed and the other with the parameters
  # can be called with a seed to reproduce a certain series of runs

  parameters = DataFrame(runs=runs, w0=wt_reproduction_rate, s_del=s_del, mu_del=mu_del, N0=founder_population_size, genmax=max_generations, SEED=rng.seed)
  ratchet_speeds = DataFrame(runs=1:runs)

  mean_speed = Vector{Float64}(undef, runs)
  for i in eachindex(run_K)
    global carrying_capacity = round(Int, 10^run_K[i]) # be careful, the carrying capacity is given by 10^(input values) rounded to an Integer
    for j = 1:runs
      mean_speed[j] = 1 / mean( deleteat!(times_until_loss_fittest_class(), 1)) # do not include the loss of the zero-mutations class
    end
    ratchet_speeds[:, "$(run_K[i])"] = mean_speed
  end
  CSV.write(string("data/ratchet_speeds_K_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), ratchet_speeds)
  CSV.write(string("data/ratchet_speeds_K_para_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), parameters)
  return ratchet_speeds, parameters
end

function ratchet_speeds_w0(run_w0::Vector{Float64}; rng=Random.seed!(), runs::Int=1)
  # runs the simulation for different wild-type reproduction rates given by 10^run_w0
  # creates two output files which can be uniquely identified by the seed, one file with the mean ratchet speed and the other with the parameters
  # can be called with a seed to reproduce a certain series of runs

  parameters = DataFrame(runs=runs, s_del=s_del, mu_del=mu_del, N0=founder_population_size, K=carrying_capacity, genmax=max_generations, SEED=rng.seed)
  ratchet_speeds = DataFrame(runs=1:runs)

  mean_speed = Vector{Float64}(undef, runs)
  for i in eachindex(run_w0)
    global wt_reproduction_rate = 10^run_w0[i] # be careful, the wild-type reporduction rate is given by 10^(input values)
    for j = 1:runs
      mean_speed[j] = 1 / mean( deleteat!(times_until_loss_fittest_class(), 1)) # do not include the loss of the zero-mutations class
    end
    ratchet_speeds[:, "$(run_w0[i])"] = mean_speed
  end
  CSV.write(string("data/ratchet_speeds_w0_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), ratchet_speeds)
  CSV.write(string("data/ratchet_speeds_w0_para_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), parameters)
  return ratchet_speeds, parameters
end

function ratchet_speeds_report()
  # the ratchet speed data used in the report was generated with this function

  run_s_del=collect(-3.:0.1:-1.5)
  run_mu_del=collect(-1.:0.05:-0.)
  init_quick(1000)
  global mu_del = 10^run_mu_del[8]
  ratchet_speeds_s_del(run_s_del, runs=50)
  init_quick(1000)
  global s_del = -10^run_s_del[8]
  ratchet_speeds_mu_del(run_mu_del, runs=50)
  init_quick()
  run_K = collect(2.:0.2:4.)
  global s_del = -0.005
  global max_generations = 20000
  ratchet_speeds_K(run_K, runs=50)
  run_w0 = collect(0.1:0.05:0.6)
  init_quick(1000)
  global s_del = -0.005
  ratchet_speeds_w0(run_w0, runs=50)

  nothing
end

function times_phases_range_s_del_mu_del(run_s_del::Vector{Float64}, run_mu_del::Vector{Float64}; rng=Random.seed!(), runs::Int=2)
  # runs the simulation for different deleterious selection coefficients given by -10^run_s_del
  # and different mutation rates given by 10^run_mu_del going from high to low mutation rates (this is important!)
  # creates seven output files which can be uniquely identified by the seed
  # for each of the three phases (equilibration, ratchet, meltdown) the mean duration and the variance
  # plus one file with the parameters
  # can be called with a seed to reproduce a certain series of runs

  parameters = DataFrame(runs=runs, w0=wt_reproduction_rate, N0=founder_population_size, K=carrying_capacity, genmax=max_generations, SEED=rng.seed)
  equilibration_time_mean = DataFrame(range=run_mu_del)
  equilibration_time_var = DataFrame(range=run_mu_del)
  ratchet_time_mean = DataFrame(range=run_mu_del)
  ratchet_time_var = DataFrame(range=run_mu_del)
  meltdown_time_mean = DataFrame(range=run_mu_del)
  meltdown_time_var = DataFrame(range=run_mu_del)

  equilibration_time_runs = Vector{Int}(undef, runs)
  ratchet_time_runs = Vector{Int}(undef, runs)
  meltdown_time_runs = Vector{Int}(undef, runs)
  survival = 0

  equilibration_time_mean_mu_del = Vector{Float64}(undef, length(run_mu_del))
  ratchet_time_mean_mu_del = Vector{Float64}(undef, length(run_mu_del))
  meltdown_time_mean_mu_del = Vector{Float64}(undef, length(run_mu_del))
  equilibration_time_var_mu_del = Vector{Float64}(undef, length(run_mu_del))
  ratchet_time_var_mu_del = Vector{Float64}(undef, length(run_mu_del))
  meltdown_time_var_mu_del = Vector{Float64}(undef, length(run_mu_del))

  for i in eachindex(run_s_del)
    global s_del = -10^run_s_del[i] # be careful, the selection coefficient is given by -10^(input values)
    for j in eachindex(run_mu_del)
      global mu_del = 10^run_mu_del[j] # be careful, the mutation rate is given by 10^(input values)
      for k = 1:runs
        survival = 0
        equilibration_time_runs[k], ratchet_time_runs[k], meltdown_time_runs[k] = times_phases()
        if equilibration_time_runs[k] + ratchet_time_runs[k] + meltdown_time_runs[k] > max_generations
          survival = 1
          break
        end
      end
      if survival > 0 # if at least one run did not result in extinction, lower mutations rates are not considered
        for k = j:length(run_mu_del)
          equilibration_time_mean_mu_del[k] = Inf
          ratchet_time_mean_mu_del[k] = Inf
          meltdown_time_mean_mu_del[k] = Inf
          equilibration_time_var_mu_del[k] = Inf
          ratchet_time_var_mu_del[k] = Inf
          meltdown_time_var_mu_del[k] = Inf
        end
        println("Simulation aborted at s_del=", run_s_del[i], " mu_del=", run_mu_del[j], ", increase max_generations!")
        break
      end
      equilibration_time_mean_mu_del[j] = mean(equilibration_time_runs)
      ratchet_time_mean_mu_del[j] = mean(ratchet_time_runs)
      meltdown_time_mean_mu_del[j] = mean(meltdown_time_runs)
      equilibration_time_var_mu_del[j] = var(equilibration_time_runs)
      ratchet_time_var_mu_del[j] = var(ratchet_time_runs)
      meltdown_time_var_mu_del[j] = var(meltdown_time_runs)
    end
    equilibration_time_mean[:, "$(run_s_del[i])"] = equilibration_time_mean_mu_del
    ratchet_time_mean[:, "$(run_s_del[i])"] = ratchet_time_mean_mu_del
    meltdown_time_mean[:, "$(run_s_del[i])"] = meltdown_time_mean_mu_del
    equilibration_time_var[:, "$(run_s_del[i])"] = equilibration_time_var_mu_del
    ratchet_time_var[:, "$(run_s_del[i])"] = ratchet_time_var_mu_del
    meltdown_time_var[:, "$(run_s_del[i])"] = meltdown_time_var_mu_del
  end

  CSV.write(string("data/equilibration_time_mean_K", carrying_capacity, "_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), equilibration_time_mean)
  CSV.write(string("data/ratchet_time_mean_K", carrying_capacity, "_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), ratchet_time_mean)
  CSV.write(string("data/meltdown_time_mean_K", carrying_capacity, "_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), meltdown_time_mean)
  CSV.write(string("data/equilibration_time_var_K", carrying_capacity, "_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), equilibration_time_var)
  CSV.write(string("data/ratchet_time_var_K", carrying_capacity, "_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), ratchet_time_var)
  CSV.write(string("data/meltdown_time_var_K", carrying_capacity, "_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), meltdown_time_var)
  CSV.write(string("data/times_phases_parameters_K", carrying_capacity, "_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), parameters)

  return equilibration_time_mean, ratchet_time_mean, meltdown_time_mean, equilibration_time_var, ratchet_time_var, meltdown_time_var, parameters
end

function times_phases_range_K(run_K::Vector{Float64}; rng=Random.seed!(), runs::Int=2)
  # runs the simulation for different carrying capacities given by 10^run_K
  # creates seven output files which can be uniquely identified by the seed
  # for each of the three phases (equilibration, ratchet, meltdown) the mean duration and the variance
  # plus one file with the parameters
  # can be called with a seed to reproduce a certain series of runs

  parameters = DataFrame(runs=runs, w0=wt_reproduction_rate, s_del=s_del, mu_del=mu_del, N0=founder_population_size, genmax=max_generations, SEED=rng.seed)
  times_mean_var = DataFrame(times=["eq_mean", "eq_var", "r_mean", "r_var", "m_mean", "m_var"])

  equilibration_time = Vector{Int}(undef, runs)
  ratchet_time = Vector{Int}(undef, runs)
  meltdown_time = Vector{Int}(undef, runs)
  survival = 0

  for i in eachindex(run_K)
    global carrying_capacity = round(Int, 10^run_K[i]) # be careful, the carrying capacity is given by 10^(input values) rounded to an Integer
    for k = 1:runs
      survival = 0
      equilibration_time[k], ratchet_time[k], meltdown_time[k] = times_phases()
      if equilibration_time[k] + ratchet_time[k] + meltdown_time[k] > max_generations
        survival = 1
        break
      end
    end
    if survival > 0 # if at least one run did not result in extinction, higher carrying capacities are not considered
      for k = i:length(run_K)
        times_mean_var[:, "$(run_K[k])"] = fill(Inf, 6)
      end
      println("Simulation aborted at K=", run_K[i], ", increase max_generations!")
      break
    end
    times_mean_var[:, "$(run_K[i])"] = [mean(equilibration_time), var(equilibration_time), mean(ratchet_time), var(ratchet_time), mean(meltdown_time), var(meltdown_time)]
  end

  CSV.write(string("data/times_phases_K_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), times_mean_var)
  CSV.write(string("data/times_phases_parameters_K_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), parameters)

  return times_mean_var, parameters
end

function times_phases_range_w0(run_w0::Vector{Float64}; rng=Random.seed!(), runs::Int=2)
  # runs the simulation for different wild-type reproduction rates given by 10^run_w0
  # creates seven output files which can be uniquely identified by the seed
  # for each of the three phases (equilibration, ratchet, meltdown) the mean duration and the variance
  # plus one file with the parameters
  # can be called with a seed to reproduce a certain series of runs

  parameters = DataFrame(runs=runs, s_del=s_del, mu_del=mu_del, N0=founder_population_size, K=carrying_capacity, genmax=max_generations, SEED=rng.seed)
  times_mean_var = DataFrame(times=["eq_mean", "eq_var", "r_mean", "r_var", "m_mean", "m_var"])

  equilibration_time = Vector{Int}(undef, runs)
  ratchet_time = Vector{Int}(undef, runs)
  meltdown_time = Vector{Int}(undef, runs)
  survival = 0

  for i in eachindex(run_w0)
    global wt_reproduction_rate = 10^run_w0[i] # be careful, the wild-type reporduction rate is given by 10^(input values)
    for k = 1:runs
      survival = 0
      equilibration_time[k], ratchet_time[k], meltdown_time[k] = times_phases()
      if equilibration_time[k] + ratchet_time[k] + meltdown_time[k] > max_generations
        survival = 1
        break
      end
    end
    if survival > 0 # if at least one run did not result in extinction, higher carrying capacities are not considered
      for k = i:length(run_w0)
        times_mean_var[:, "$(run_w0[k])"] = fill(Inf, 6)
      end
      println("Simulation aborted at w0=", run_w0[i], ", increase max_generations!")
      break
    end
    times_mean_var[:, "$(run_w0[i])"] = [mean(equilibration_time), var(equilibration_time), mean(ratchet_time), var(ratchet_time), mean(meltdown_time), var(meltdown_time)]
  end

  CSV.write(string("data/times_phases_w0_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), times_mean_var)
  CSV.write(string("data/times_phases_parameters_w0_", parameters.SEED[1], parameters.SEED[2], parameters.SEED[3], parameters.SEED[4]), parameters)

  return times_mean_var, parameters
end

function times_phases_report()
   # the data used in the report was generated with this function

  run_s_del=collect(-3.:0.1:-1.5)
  run_mu_del=collect(-0.:-0.05:-1.) # going from high to low mutation rates (this is important for plotting!)
  init_quick(100)
  global max_generations = 20000
  times_phases_range_s_del_mu_del(run_s_del, run_mu_del, runs=150)
  init_quick(1000)
  global max_generations = 20000
  times_phases_range_s_del_mu_del(run_s_del, run_mu_del, runs=50)
  init_quick(10000)
  global max_generations = 20000
  times_phases_range_s_del_mu_del(run_s_del, run_mu_del, runs=20)
  run_K = collect(2.:0.2:4.)
  init_quick()
  global max_generations = 20000
  global s_del = -0.005
  times_phases_range_K(run_K, runs=50)
  run_w0 = collect(0.1:0.05:0.6)
  init_quick(1000)
  global max_generations = 20000
  global s_del = -0.005
  times_phases_range_w0(run_w0, runs=50)

  nothing
end
