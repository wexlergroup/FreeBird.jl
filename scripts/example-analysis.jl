df = read_output("scripts/LJ_111_p08_output_df.csv")

gi = ωᵢ(df.iter, 640)

ei = df.emax .- minimum(df.emax)
ts = collect(10:1:2000)
kb = 8.617333262e-5 # eV/K
eps = 0.1
beta = 1 ./(kb.*ts)
dof = 24

zs = [partition_function(b, gi, ei) for b in beta]

u = [internal_energy(b, gi, ei) for b in beta]

# cvs = [cv(b, gi, ei, dof) for b in beta]

cvs = cv(df, beta, dof, 640)

using Plots

plot(ts, cvs, label="Cv-FreeBird", xlabel="T(K)", ylabel="Heat Capacity", title="LJ(111) 4x4 θ=8/16")