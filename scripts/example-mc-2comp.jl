# Define parameters for MCMC simulation
L = 4 #s quare_supercell_dimensions[1]
M = 4 # square_supercell_dimensions[2]
Ns = [4]
adsorption_energy = -0.04
nn_energy = -0.01
nnn_energy = -0.0025
temperatures = 200:-1:1
num_equilibration_steps = 10_000
num_sampling_steps = 10_000
random_seed = 1234
k_B = 8.617_333_262e-5


df_samples = DataFrame(N = Int[], T = Float64[], E = Float64[], Cv = Float64[], acceptance_rate = Float64[])

# Perform the Monte Carlo simulation
for N in Ns
    println("N = $N -------------------")

    # Initialize the lattice
    square_occupations = [false for i in 1:square_supercell_dimensions[1]*square_supercell_dimensions[2]*length(square_basis)]
    for i in sample(1:length(square_occupations), N, replace=false)
        square_occupations[i] = true
    end

    initial_lattice = MLattice{2,SquareLattice}(;
        basis=square_basis,
        supercell_dimensions = square_supercell_dimensions,
    )
    initial_lattice.components[1][1:6] .= false
    initial_lattice.components[2][11:16] .= false

    h = MLatticeHamiltonian(2,[GenericLatticeHamiltonian(adsorption_energy, [nn_energy*i^2, nnn_energy], u"eV") for i in 1:3])

    for temp in temperatures
        println("T = $temp")

        # Equilibrate the lattice
        temperature = Float64(temp)
        equilibration_energies, equilibration_configurations, equilibration_accepted_steps = nvt_monte_carlo(
            initial_lattice,
            h,
            temperature,
            num_equilibration_steps,
            random_seed
        )

        # Sample the lattice
        sampling_energies, sampling_configurations, sampling_accepted_steps = nvt_monte_carlo(
            equilibration_configurations[end],
            h,
            temperature,
            num_sampling_steps,
            random_seed
        )

        # Compute the heat capacity
        E = mean(sampling_energies)
        Cv = var(sampling_energies) / (k_B * temperature^2)

        # Compute the acceptance rate
        acceptance_rate = sampling_accepted_steps / num_sampling_steps

        # Append the results to the DataFrame
        append!(df_samples, DataFrame(N = N, T = temperature, E = E, Cv = Cv, acceptance_rate = acceptance_rate))

        # Update the initial lattice
        initial_lattice = deepcopy(sampling_configurations[end])
    end
end