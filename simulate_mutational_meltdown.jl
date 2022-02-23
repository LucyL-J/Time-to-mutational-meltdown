import Pkg
Pkg.activate("packages")
Pkg.instantiate()
using Distributions, StatsBase, Random, DataFrames, CSV

function init_quick()

  global founder_population_size  = 20
  global carrying_capacity        = 10^3
  global max_generations          = 1000 # if the population did not get extinct after this number of generations, the simulation is stopped
  global wt_reproduction_rate     = 2.0
  global mu                       = 0.27 # mutation rate per genome per generation
  global s                        = -0.005 # selection coeffiecient is negative here
  global random_seed              = abs(rand(Int)) # for random number generation

  nothing
end

function init_quick(K::Int) # to initialise with the standard parameters but a different carrying capacity

  init_quick()
  global carrying_capacity = K

  nothing
end

function create_founder_population()

  population = zeros(Int, founder_population_size)

end

function mutation!(population)

  for indiv_id in eachindex(population)

    new_del_mut = rand( Poisson( mu ) )

    # if new mutations have occured, update individual
    if new_del_mut != 0
      population[indiv_id] += new_del_mut
    end
  end

  nothing
end

function clonal_proliferation!(population)

  list_of_indiv_without_offspring = Array{Int}(undef, 0)

  for indiv_id in eachindex(population)
    num_offspring = rand( Poisson(  wt_reproduction_rate * ( 1 + s )^population[indiv_id] ) )
    if num_offspring < 1
      push!( list_of_indiv_without_offspring, indiv_id )
    elseif num_offspring > 1
      for n = 1:num_offspring-1
      push!( population, population[indiv_id] )
      end
    end
  end
  deleteat!( population, list_of_indiv_without_offspring ) # individuals which do not have offspring are removed from the population

  nothing
end

function apply_carrying_capacity!(population)

  if ( ( length(population) - carrying_capacity ) > 0 ) # effect only if the population size actually exceeds the carrying capacity
    # individuals are deleted uniformly at random
    deleteat!( population, sort!(sample( 1:length(population), length(population) - carrying_capacity, replace=false) ) )
  end

  nothing
end

function extinction_time() # returns the extinction time

  population = create_founder_population( )

  @inbounds for generation in 1:max_generations
    mutation!(population)
    clonal_proliferation!(population)
    apply_carrying_capacity!(population)

    if length(population) == 0
      return generation
      break
    end
  end

  return max_generations+1 # if the population did not get extinct
end

function num_del_mutations(population)

    count_mutations = Vector{Int}(undef, 0)
    for indiv_id in eachindex(population)
        if population[indiv_id]+1 > length(count_mutations)
            for i = 1:population[indiv_id]+1 - length(count_mutations)
                push!(count_mutations, 0)
            end
        end
        count_mutations[population[indiv_id]+1] += 1
    end

    count_mutations
end

function add_gen!(distribution_mutations, dist_generation, generation)

    if length(dist_generation) > size(distribution_mutations, 1)
        for i = size(distribution_mutations, 1):length(dist_generation)-1
            push!(distribution_mutations, [i; zeros(Int, size(distribution_mutations, 2)-1)])
        end
    else
        for i = length(dist_generation):size(distribution_mutations, 1)-1
            push!(dist_generation, 0)
        end
    end

    distribution_mutations[!, "gen$generation"] = dist_generation
end

function extinction_time_parameters(;rng=Random.seed!(random_seed))
  # returns a DataFrame with the paramters used for the simulation (including the seed) plus the extinction time T
  # can be called with a seed to reproduce a certain run

  parameters = DataFrame(w0=wt_reproduction_rate, s=s, mu=mu, N0=founder_population_size, K=carrying_capacity, genmax=max_generations, seed=random_seed)

  population = create_founder_population( )

  @inbounds for generation in 1:max_generations
    mutation!(population)
    clonal_proliferation!(population)
    apply_carrying_capacity!(population)

    if length(population) == 0
      return parameters
      break
    end
  end

  return parameters # has four rows only because the of the size of the seed
end

function extinction_time_dist_mutations(;rng=Random.seed!(random_seed))
  # returns two DataFrames
  # parameters with the paramters used for the simulation (including the raondom seed) plus the extinction time T
  # distribution_mutations: number of individuals with certain number of mutations (row) at a certain generation (column)
  # can be called with a random seed to reproduce a certain run

  parameters = DataFrame(w0=wt_reproduction_rate, s=s, mu=mu, N0=founder_population_size, K=carrying_capacity, genmax=max_generations, seed=random_seed)

  population = create_founder_population( )

  distribution_mutations = DataFrame(del_muts = [0], gen0 = [founder_population_size])

  @inbounds for generation in 1:max_generations
    mutation!(population)
    clonal_proliferation!(population)
    apply_carrying_capacity!(population)

    if length(population) == 0
      return distribution_mutations, parameters
      break
    end

    add_gen!(distribution_mutations, num_del_mutations(population), generation)
  end

  return distribution_mutations, parameters
end

function times_until_loss_fittest_class()

    critical_num_mutations = round(Int, -log(wt_reproduction_rate)/s, RoundUp) # stop when the critical number of mutations is reached and the actual meltdown starts
    tau = Array{Int}(undef, 0)

    population = create_founder_population( )

    fittest_class = 0
    t = 0

    for generation in 1:max_generations
      mutation!(population)
      clonal_proliferation!(population)
      apply_carrying_capacity!(population)

      if minimum(population) > fittest_class
        append!(tau, generation - t)
        t = generation
        fittest_class = minimum(population)
      end

      if mean(population) >= critical_num_mutations
        break
      end

    end

    if length(tau) < 1 # the zero-mutation class is never lost
      return [Inf]
    else
      return tau
    end
  end

  function times_phases()

    critical_num_mutations = round(Int, -log(wt_reproduction_rate)/s, RoundUp)
    equilibration_time, ratchet_time = 0, 0

    population = create_founder_population( )

    for generation in 1:max_generations
      mutation!(population)
      clonal_proliferation!(population)
      apply_carrying_capacity!(population)

      if length(population) == 0
        meltdown_time = generation - equilibration_time - ratchet_time
        return equilibration_time, ratchet_time, meltdown_time
        break
      elseif minimum(population) <= 0
        equilibration_time = generation
      elseif mean(population) <= critical_num_mutations
        ratchet_time = generation - equilibration_time
      end
    end

    return equilibration_time, ratchet_time, max_generations+1-equilibration_time-ratchet_time
  end
