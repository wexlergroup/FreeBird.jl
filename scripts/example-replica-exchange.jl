# Initialize replicas with random configurations
num_replicas = 5
replicas = [deepcopy(initial_lattice) for _ in 1:num_replicas]
for replica in replicas
    replica.occupations = [false for i in 1:length(replica.occupations)]
    for i in sample(1:length(replica.occupations), 4, replace=false)
        replica.occupations[i] = true
    end
end

# Define parameters for replica exchange simulation
temperatures = [50., 100., 150., 200., 250.]  # K
num_steps = 100  # 10000
swap_fraction = 0.1  # 10% swap, 90% displacement
k_B = 8.617_333_262e-5
for replica in replicas
    println(replica.occupations)
end

# Run the simulation
replicas, energies, temperature_history = SamplingSchemes.nvt_replica_exchange(
    replicas,
    temperatures,
    num_steps,
    ham,
    swap_fraction
)

# Plot the temperature history
p = plot(title="Temperature History", xlabel="Step", ylabel="Temperature", legend=:outerright)
for i in 1:num_replicas
    plot!(p, temperature_history[:, i], label="Replica $i")
end
display(p)