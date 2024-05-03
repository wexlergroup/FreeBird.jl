# First, read in the output energies and iterations as a DataFrame
df = read_output("scripts/LJ_111_p08_output_df.csv")

# Calculate the ω factors
ωi = ωᵢ(df.iter, 640)

# Shift the energies to be greater than or equal to zero
Ei = df.emax .- minimum(df.emax)
# Specify the temperatures that we are interested in
Ts = collect(10:1:2000)
# Define the Boltzmann constant
kb = 8.617333262e-5 # eV/K
# Calculate the inverse temperatures
β = 1 ./(kb.*Ts)
# Define the degrees of freedom, which is 3×8 for the 8-particle system
dof = 24

# Calculate the partition functions for each temperature
Zs = [partition_function(b, ωi, Ei) for b in β]

# Calculate the internal energies for each temperature
U = [internal_energy(b, ωi, Ei) for b in β]

# Calculate the heat capacities as a function of temperature
cvs = cv(df, β, dof, 640)

# Plot the heat capacities
using Plots
plot(Ts, cvs, label="\$C_V\$", xlabel="\$T(K)\$", ylabel="Heat Capacity", title="LJ(111) 4x4 θ=8/16")