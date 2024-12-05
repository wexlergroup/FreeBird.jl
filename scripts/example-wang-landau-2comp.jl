using FreeBird, Distributions, DataFrames, Plots

square_supercell_dimensions = (4, 4, 1)
square_basis=[(0.0, 0.0, 0.0)]
square_occupations = [false for i in 1:square_supercell_dimensions[1]*square_supercell_dimensions[2]*length(square_basis)]

adsorption_energy = -0.04
nn_energy = -0.01
nnn_energy = -0.0025

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

# Initialize the lattice
occupied_sites = sample(1:length(square_occupations), 4, replace=false)

# initial_lattice = SLattice{SquareLattice}(;
#            supercell_dimensions = square_supercell_dimensions,
#            occupations=occupied_sites)

initial_lattice = MLattice{2,SquareLattice}()
initial_lattice.components[1][1:6] .= false
initial_lattice.components[2][11:16] .= false

h = MLatticeHamiltonian(2,[GenericLatticeHamiltonian(adsorption_energy, [nn_energy*i^2, nnn_energy], u"eV") for i in 1:3])

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

# Plot g(E)
p1 = plot(E, g, xlabel="Energy (eV)", ylabel="Density of states", legend=false, title="Shifted entropy")

# Plot exp.(df_entropy.entropy)
p2 = plot(E, exp.(df_wl.entropy), xlabel="Energy (eV)", ylabel="Density of states", legend=false, title="Unshifted entropy")

# Display the plots
plot(p1, p2, layout=(1, 2), size=(800, 400))