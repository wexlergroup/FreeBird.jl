using FreeBird, Distributions, DataFrames, Plots

square_supercell_dimensions = (4, 4, 1)
square_basis=[(0.0, 0.0, 0.0)]

adsorption_energy = -0.04
nn_energy = -0.01
nnn_energy = -0.0025

# Initialize the lattice
occupied_sites = sample(1:16, 4, replace=false)

initial_lattice = LatticeSystem{SquareLattice}(;
           supercell_dimensions = square_supercell_dimensions,
           occupations=occupied_sites)

h = GenericLatticeHamiltonian(adsorption_energy, [nn_energy, nnn_energy], u"eV")


walkers = [deepcopy(initial_lattice) for i in 1:2000]

for walker in walkers
    walker.occupations = [false for i in 1:length(walker.occupations)]
    for i in sample(1:length(walker.occupations), 4, replace=false)
        walker.occupations[i] = true
    end
end

unique!(x -> x.occupations, walkers) # remove duplicates

# h = LatticeGasHamiltonian([-0.04,-0.01,-0.0025].*u"eV"...)
# h = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

# Now, we can create the liveset by combining the walkers and the potential
ls = LatticeGasWalkers(deepcopy(LatticeWalker.(walkers[1:1000])), h, perturb_energy=1e-10)

# We can now create a Monte Carlo object
mc = MCNewSample()

# Let's specify the nested sampling parameters: 
# Monte Carlo steps, energy perturbation, current fail count, allowed fail count
# Use `?LatticeNestedSamplingParameters` to see the documentation for more information
ns_params = LatticeNestedSamplingParameters(1,1e-10,0,1000000)
# We can also create a save object, using default value for most parameters except `n_snap` for saving snapshots
save = SaveEveryN(n_traj=10000, n_snap=10000)

# And finally, we can run the simulation
ns_energies, ls, _ = nested_sampling_loop!(ls, ns_params, 500000, mc, save; accept_same_config=true)

num_walkers = 1000
Ts = collect(1.0:0.1:200.0)
k_B = 8.617_333_262e-5

# Compute ω_i
ω_0 = 16 * 15 * 14 * 13 / 4 / 3 / 2 / 1  # 16 choose 4
ω_i = ω_0 * (1 / num_walkers) * (num_walkers / (num_walkers + 1)) .^ (collect(1:length(ns_energies.emax)) .- 1)

# Compute the partition function as a function of temperature
Z = zeros(length(Ts))
E_i_ns = ns_energies.emax
E_rel_i_ns = E_i_ns .- minimum(E_i_ns)

for (i, temp) in enumerate(Ts)
    Z[i] = sum(exp.(-E_rel_i_ns ./ (k_B * temp)) .* ω_i)
end

# Compute the average energy as a function of temperature
E_avg_i_ns = zeros(length(Ts))
for (i, temp) in enumerate(Ts)
    E_avg_i_ns[i] = sum(E_rel_i_ns .* exp.(-E_rel_i_ns ./ (k_B * temp)) .* ω_i) / Z[i]
end

# Compute the average energy squared as a function of temperature
E2_avg_i_ns = zeros(length(Ts))
for (i, temp) in enumerate(Ts)
    E2_avg_i_ns[i] = sum(E_rel_i_ns.^2 .* exp.(-E_rel_i_ns ./ (k_B * temp)) .* ω_i) / Z[i]
end

# Compute the heat capacity as a function of temperature
Cv_i_ns = (E2_avg_i_ns .- E_avg_i_ns.^2) ./ (k_B * Ts.^2) / k_B

# Plot the heat capacity as a function of temperature
p = plot(Ts, Cv_i_ns, xlabel="Temperature (K)", ylabel="Heat Capacity (k_B)", label="NS", title="Heat Capacity vs. Temperature")

# Wang-Landau

# Define parameters for Wang-Landau simulation
num_steps = 100
flatness_criterion = 0.8
f_initial = Float64(MathConstants.e)
f_min = exp(10^-8)
energy_min = 20.5 * nn_energy + nn_energy / 8
energy_max = 16 * nn_energy - nn_energy / 8
num_energy_bins = 100
energy_bins = collect(range(energy_min, stop=energy_max, length=num_energy_bins))
random_seed = 1234


entropy, histogram, bin_energies, energies, configurations = wang_landau(
                    initial_lattice,
                    h,
                    num_steps,
                    flatness_criterion,
                    f_initial,
                    f_min,
                    energy_bins,
                    random_seed
)

df_entropy = DataFrame(energy = bin_energies, entropy = entropy)

df_wl = df_entropy
S_shifted = df_wl.entropy .- minimum(df_wl.entropy[df_wl.entropy .> 0])
g = exp.(S_shifted)
E = df_wl.energy



# Compute the partition function as a function of temperature
E_rel = E .- minimum(E)
Ts = collect(1.0:0.1:200.0)
Z = zeros(length(Ts))
for (i, temp) in enumerate(Ts)
    Z[i] = sum(exp.(-E_rel ./ (k_B * temp)) .* g)
end

# Plot the partition function
# plot(Ts, Z, xlabel="Temperature (K)", ylabel="Partition function", legend=false, title="Partition function")

# Compute the average energy as a function of temperature
E_avg = zeros(length(Ts))
for (i, temp) in enumerate(Ts)
    E_avg[i] = sum(E_rel .* exp.(-E_rel ./ (k_B * temp)) .* g) / Z[i]
end

# Plot the average energy
# plot(Ts, E_avg, xlabel="Temperature (K)", ylabel="Average energy (eV)", legend=false, title="Average energy")

# Compute the average energy squared as a function of temperature
E2_avg = zeros(length(Ts))
for (i, temp) in enumerate(Ts)
    E2_avg[i] = sum(E_rel.^2 .* exp.(-E_rel ./ (k_B * temp)) .* g) / Z[i]
end

# Compute the heat capacity as a function of temperature
Cv = (E2_avg .- E_avg.^2) ./ (k_B * Ts.^2) / k_B

# Plot the heat capacity as a function of temperature
plot!(p, Ts, Cv, xlabel="Temperature (K)", ylabel="Heat Capacity (k_B)", label="WL", title="Heat Capacity vs. Temperature")


# Exact enumeration

cutoff_radii = [1.1, 1.5]  # Angstrom

energies, configurations, walkers = exact_enumeration(
                          initial_lattice,
                          cutoff_radii,
                          h)

                                # Define temperature range and corresponding beta values
Ts = collect(1.0:0.1:200.0)  # Temperatures in K
βs = 1.0 ./ (k_B * Ts)  # 1/eV

# Prepare energy values for Cv calculation
ωis = ones(length(energies))
Eis = [ustrip(energy) for energy in energies]
dof = 0

# Calculate Cv for each temperature
Cv_values = [cv(β, ωis, Eis, dof) / k_B for β in βs]  # Cv in units of kB

# Plot the heat capacity as a function of temperature
plot!(p, Ts, Cv_values, xlabel="Temperature (K)", ylabel="Heat Capacity (k_B)", label="Exact", title="Heat Capacity vs. Temperature")