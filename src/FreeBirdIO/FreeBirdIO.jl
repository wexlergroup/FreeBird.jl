"""
    FreeBirdIO

Module for input/output operations in the FreeBird package.
"""
module FreeBirdIO

using AtomsIO, ExtXYZ, DataFrames, CSV
using AtomsBase, Unitful
using Arrow
using ASEconvert, PythonCall

using ..AbstractPotentials
using ..AbstractWalkers
using ..AbstractLiveSets
using ..EnergyEval

export read_single_config, read_configs, read_single_walker, read_walkers
export write_single_walker, write_walkers

export append_system

export DataSavingStrategy, SaveEveryN, SaveFreePartEveryN
export write_df, write_df_every_n, write_walker_every_n, write_ls_every_n

export generate_initial_configs, generate_accessible_configs
export convert_system_to_walker, convert_walker_to_system

# include the save strategies
include("save_strategies.jl")

"""
    set_pbc(at::Atoms, pbc::Vector)

Set the periodic boundary conditions for a system of atoms.

# Arguments
- `at::Atoms`: The system of atoms.
- `pbc::Vector`: A vector of length 3 specifying the periodic boundary conditions for each dimension. Each element can be either `true` for periodic boundary conditions or `false` for Dirichlet zero boundary conditions.

# Returns
- `FlexibleSystem`: A flexible system with the specified boundary conditions.

"""
function set_pbc(at::Atoms, pbc::Vector{Bool})
    return FlexibleSystem(at, periodicity=Tuple(pbc))
end

@inline function _has_atomic_lattice_metadata(data)
    return haskey(data, :atomic_lattice_schema) ||
           get(data, :freebird_walker, "") == "AtomicLattice"
end

"""
    convert_system_to_walker(at::FlexibleSystem, resume::Bool)

Converts a `FlexibleSystem` to a walker, using its representation metadata.

# Arguments
- `at::FlexibleSystem`: The `FlexibleSystem` object to convert.
- `resume::Bool`: Whether to resume from previous data.

# Returns
- `AtomWalker` or `LatticeWalker`: The converted walker. AtomicLattice
  metadata selects a `LatticeWalker`; ordinary atomic frames select an
  `AtomWalker`.

"""
function convert_system_to_walker(at::FlexibleSystem, resume::Bool)
    data = at.data
    if _has_atomic_lattice_metadata(data)
        return convert_system_to_atomic_lattice_walker(at, resume)
    end
    energy = (haskey(data, :energy) && resume) ? data[:energy]*u"eV" : 0.0u"eV"
    iter = (haskey(data, :iter) && resume && data[:iter]>=0) ? data[:iter] : 0
    list_num_par = (haskey(data, :list_num_par) && resume) ? data[:list_num_par] : [length(at)]
    C = length(list_num_par)
    frozen = (haskey(data, :frozen) && resume) ? data[:frozen] : zeros(Bool, C)
    e_frozen = (haskey(data, :energy_frozen_part) && resume) ? data[:energy_frozen_part]*u"eV" : 0.0u"eV"
    new_list = [Atom(atomic_symbol(i),position(i)) for i in at.particles]
    at = FastSystem(new_list, cell_vectors(at), periodicity(at))
    return AtomWalker{C}(at; energy=energy, iter=iter, list_num_par=list_num_par, frozen=[frozen...], energy_frozen_part=e_frozen)
end

_split_metadata(value::AbstractVector) = string.(value)
_split_metadata(value::Tuple) = string.(value)
_split_metadata(value) = isempty(string(value)) ? String[] : String.(split(string(value), ','))

# A non-numeric prefix prevents extXYZ's fixed-decimal float formatting from
# rounding the constructor scalars needed to reproduce the site fingerprint.
_encode_atomic_geometry_scalar(value::Float64) = "f64=" * repr(value)
_decode_atomic_geometry_scalar(value::Real) = Float64(value)
function _decode_atomic_geometry_scalar(value::AbstractString)
    startswith(value, "f64=") || throw(ArgumentError(
        "AtomicLattice geometry scalars must be numeric or f64-prefixed text"))
    return parse(Float64, value[5:end])
end

# ExtXYZ treats the standalone element symbol F as a Boolean token. Prefix
# chemical-symbol metadata so every element is preserved as text.
_encode_atomic_symbols(symbols) = "symbols=" * join(symbols, ',')
_decode_atomic_symbols(value::Bool) = [value ? "T" : "F"]
function _decode_atomic_symbols(value)
    if value isa AbstractString && startswith(value, "symbols=")
        return _split_metadata(value[9:end])
    end
    return _split_metadata(value)
end

function _atomic_lattice_geometry(name::AbstractString)
    name == "SquareLattice" && return SquareLattice
    name == "TriangularLattice" && return TriangularLattice
    name == "GenericLattice" && return GenericLattice
    throw(ArgumentError("unsupported AtomicLattice geometry '$name' in trajectory"))
end

# Stable, process-independent checksum of the ordered adsorption-site
# coordinates. This detects trajectory/site-order incompatibilities caused by
# a changed ASE builder without storing a second full copy of the coordinates.
function _atomic_site_geometry_fingerprint(lattice::AtomicLattice)
    state = UInt64(0xcbf29ce484222325)
    function mix_word(state::UInt64, word::UInt64)
        for shift in 0:8:56
            state = (state ⊻ ((word >> shift) & 0xff)) * UInt64(0x100000001b3)
        end
        return state
    end
    function mix_coordinate(state::UInt64, coordinate::Real)
        scaled = Float64(coordinate) * 1e10
        if isfinite(scaled) && -0x1p63 <= scaled < 0x1p63
            # Use fixed-point encoding for coordinates within the Int64 range.
            return mix_word(state, reinterpret(UInt64, round(Int64, scaled)))
        end
        # A finite coordinate can overflow the fixed-point Int64 encoding.
        # Mark that alternate representation explicitly before mixing its
        # rounded floating-point bits to keep the checksum deterministic.
        state = mix_word(state, UInt64(0x6f766572666c6f77)) # "overflow"
        rounded = round(Float64(coordinate); digits=10)
        rounded == 0.0 && (rounded = 0.0) # canonicalize signed zero
        return mix_word(state, reinterpret(UInt64, rounded))
    end
    state = mix_word(state, UInt64(length(lattice.all_sites)))
    for (x, y) in lattice.all_sites
        # Ten decimal places are far below the site's physical tolerance but
        # avoid rejecting a trajectory for harmless last-bit ASE drift.
        state = mix_coordinate(state, x)
        state = mix_coordinate(state, y)
    end
    # The non-numeric prefix forces ExtXYZ to preserve the checksum as text.
    return "fnv64=" * string(state; base=16, pad=16)
end

_angstrom_values(vector) = Float64[ustrip(u"Å", value) for value in vector]

"""Check that extXYZ atoms and cell agree with the reconstructing metadata."""
function _validate_atomic_serialized_frame(at::FlexibleSystem,
                                           lattice::AtomicLattice)
    expected = pyconvert(AbstractSystem, lattice.ase_lattice)
    length(at) == length(expected) || throw(ArgumentError(
        "AtomicLattice trajectory contains $(length(at)) atoms, but its " *
        "metadata reconstructs $(length(expected))"))
    stored_periodicity = haskey(at.data, :periodicity) ?
        Tuple(at.data[:periodicity]) : periodicity(at)
    stored_periodicity == periodicity(expected) || throw(ArgumentError(
        "AtomicLattice trajectory periodicity disagrees with its metadata"))
    for (stored_vector, expected_vector) in
        zip(cell_vectors(at), cell_vectors(expected))
        isapprox(_angstrom_values(stored_vector),
                 _angstrom_values(expected_vector); rtol=1e-10, atol=1e-8) ||
            throw(ArgumentError(
                "AtomicLattice trajectory cell disagrees with its metadata"))
    end
    for (stored_atom, expected_atom) in zip(at.particles, expected.particles)
        atomic_symbol(stored_atom) == atomic_symbol(expected_atom) ||
            throw(ArgumentError(
                "AtomicLattice trajectory atom ordering or species disagrees " *
                "with its metadata"))
        isapprox(_angstrom_values(position(stored_atom)),
                 _angstrom_values(position(expected_atom));
                 rtol=1e-10, atol=1e-8) || throw(ArgumentError(
            "AtomicLattice trajectory atom positions disagree with its metadata"))
    end
    return nothing
end

"""Reconstruct an `AtomicLattice` walker from FreeBird extXYZ metadata."""
function convert_system_to_atomic_lattice_walker(at::FlexibleSystem, resume::Bool)
    data = at.data
    required = (:freebird_walker, :atomic_lattice_schema,
                :lattice_atom, :adsorbate_atoms,
                :supercell_dimensions, :lattice_constant,
                :lattice_periodicity, :num_nearest_neighbors, :type_of_sites,
                :adsorbate_height, :lattice_geometry, :ase_surface,
                :image_multiplicity, :site_geometry_fingerprint,
                :component_occupations, :energy_convention, :energy, :iter)
    missing_keys = filter(k -> !haskey(data, k), required)
    isempty(missing_keys) || throw(ArgumentError(
        "AtomicLattice trajectory is missing metadata: $(join(missing_keys, ", "))"))
    schema = Int(data[:atomic_lattice_schema])
    schema == 1 || throw(ArgumentError(
        "unsupported AtomicLattice trajectory schema $schema; expected 1"))
    string(data[:freebird_walker]) == "AtomicLattice" || throw(ArgumentError(
        "AtomicLattice trajectory has invalid freebird_walker metadata; " *
        "expected 'AtomicLattice'"))
    convention = string(data[:energy_convention])
    convention == "bare_E_v1" || throw(ArgumentError(
        "unsupported AtomicLattice energy convention '$convention'; " *
        "expected 'bare_E_v1'"))

    adsorbates = _decode_atomic_symbols(data[:adsorbate_atoms])
    site_types = _split_metadata(data[:type_of_sites])
    dims = Tuple(parse.(Int, _split_metadata(data[:supercell_dimensions])))
    pbc = Tuple(Bool(parse(Int, x)) for x in _split_metadata(data[:lattice_periodicity]))
    length(dims) == 3 || throw(ArgumentError("supercell_dimensions must contain three integers"))
    length(pbc) == 3 || throw(ArgumentError("lattice_periodicity must contain three booleans"))
    all(>(0), dims) || throw(ArgumentError(
        "supercell_dimensions must contain positive integers"))
    # Every supported elemental builder supplies nx*ny*nz substrate atoms.
    # Reject impossible metadata before asking ASE to allocate that slab.
    substrate_count = try
        foldl(Base.Checked.checked_mul, dims; init=1)
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("trajectory supercell_dimensions overflow the substrate count"))
    end
    substrate_count <= length(at) || throw(ArgumentError(
        "trajectory dimensions require $substrate_count substrate atoms, " *
        "but the stored frame contains only $(length(at)) atoms"))

    geometry = _atomic_lattice_geometry(string(data[:lattice_geometry]))
    surface = string(data[:ase_surface])
    raw_multiplicity = data[:image_multiplicity]
    image_multiplicity = raw_multiplicity isa Bool ? raw_multiplicity :
        Bool(parse(Int, string(raw_multiplicity)))
    lattice_type = AtomicLattice{length(adsorbates),geometry}
    lattice_constant_c = haskey(data, :lattice_constant_c) ?
        _decode_atomic_geometry_scalar(data[:lattice_constant_c]) : nothing
    lattice = lattice_type(
        lattice_atom=only(_decode_atomic_symbols(data[:lattice_atom])),
        surface=surface,
        supercell_dimensions=dims,
        lattice_constant=_decode_atomic_geometry_scalar(data[:lattice_constant]),
        lattice_constant_c=lattice_constant_c,
        periodicity=pbc,
        adsorbate_atoms=adsorbates,
        coverage=0.0,
        num_nearest_neighbors=Int(data[:num_nearest_neighbors]),
        image_multiplicity=image_multiplicity,
        type_of_sites=site_types,
        adsorbate_height=_decode_atomic_geometry_scalar(data[:adsorbate_height]))
    expected = string(data[:site_geometry_fingerprint])
    actual = _atomic_site_geometry_fingerprint(lattice)
    expected == actual || throw(ArgumentError(
        "AtomicLattice trajectory adsorption-site geometry does not match " *
        "the sites generated by the current ASE version"))

    encoded_components = split(string(data[:component_occupations]), ';')
    length(encoded_components) == length(adsorbates) || throw(ArgumentError(
        "trajectory contains $(length(encoded_components)) component masks for " *
        "$(length(adsorbates)) adsorbate species"))
    for (component, raw) in zip(lattice.components, encoded_components)
        encoded = string(raw)
        startswith(encoded, "bits=") && (encoded = encoded[6:end])
        all(c -> c == '0' || c == '1', encoded) || throw(ArgumentError(
            "component occupation masks must contain only zeroes and ones"))
        mask = [c == '1' for c in encoded]
        length(mask) == num_sites(lattice) || throw(ArgumentError(
            "occupation mask has $(length(mask)) sites, expected $(num_sites(lattice))"))
        component .= mask
    end
    any(sum(component[site] for component in lattice.components) > 1
        for site in 1:num_sites(lattice)) && throw(ArgumentError(
            "AtomicLattice component masks overlap at one or more sites"))
    lattice.ase_dirty = true
    sync_ase_lattice!(lattice)
    _validate_atomic_serialized_frame(at, lattice)

    energy = (resume && haskey(data, :energy)) ? Float64(data[:energy]) * u"eV" : 0.0u"eV"
    iter = (resume && haskey(data, :iter) && data[:iter] >= 0) ? Int(data[:iter]) : 0
    return LatticeWalker(lattice; energy=energy, iter=iter)
end

"""
    read_single_config(filename::String, pbc::Vector)

Reads a single configuration from the specified file and sets the periodic
boundary conditions (PBC) for the atoms. Serialized `AtomicLattice` frames keep
their stored PBC so their slab geometry can be reconstructed exactly.

# Arguments
- `filename::String`: The name of the file to read the configuration from.
- `pbc::Vector`: A vector specifying the periodic boundary conditions for
  ordinary atomic configurations. It is ignored for serialized `AtomicLattice`
  frames.

# Returns
- `at::Atoms`: The atoms with the PBC set.

"""
function read_single_config(filename::String, pbc::Vector)
    at = Atoms(read_frame(filename::String))
    frame_pbc = _has_atomic_lattice_metadata(at.system_data) ?
        collect(periodicity(at)) : pbc
    return set_pbc(at, frame_pbc)
end

"""
    read_single_config(filename::String; pbc::String="TTT")

Reads a single configuration from the specified file.

# Arguments
- `filename::String`: The name of the file to read from.
- `pbc::String="TTT"`: The periodic boundary conditions. Default is "TTT".

# Returns
- The configuration read from the file.

"""
function read_single_config(filename::String; pbc::String="TTT")
    pbc = [pbc[1] == 'T', pbc[2] == 'T', pbc[3] == 'T']
    read_single_config(filename, pbc)
end

"""
    read_configs(filename::String, pbc::Vector)

Reads atomic configurations from a file and applies periodic boundary
conditions. Serialized `AtomicLattice` frames keep their stored PBC so their
slab geometry can be reconstructed exactly.

# Arguments
- `filename::String`: The name of the file containing the atomic configurations.
- `pbc::Vector`: A vector specifying the periodic boundary conditions for
  ordinary atomic configurations. It is ignored for serialized `AtomicLattice`
  frames.

# Returns
An array of atomic configurations with periodic boundary conditions applied,
except that serialized `AtomicLattice` frames retain their stored PBC.

"""
function read_configs(filename::String, pbc::Vector)
    ats = Atoms.(read_frames(filename::String))
    return [set_pbc(at,
        _has_atomic_lattice_metadata(at.system_data) ?
            collect(periodicity(at)) : pbc) for at in ats]
end

"""
    read_configs(filename::String; pbc::String="TTT")

Reads configurations from a file.

# Arguments
- `filename::String`: The name of the file to read configurations from.
- `pbc::String="TTT"`: Periodic boundary conditions. A string of length 3, where each character represents whether the corresponding dimension has periodic boundary conditions ('T') or not ('F').

# Returns
- The configurations read from the file.

"""
function read_configs(filename::String; pbc::String="TTT")
    pbc = [pbc[1] == 'T', pbc[2] == 'T', pbc[3] == 'T']
    read_configs(filename, pbc)
end

"""
    read_single_walker(filename::String; pbc::String="TTT", resume::Bool=true)

Reads a single walker from the specified file.

# Arguments
- `filename::String`: The path to the file containing the walker data.
- `pbc::String`: (optional) The periodic boundary conditions. Default is "TTT".
- `resume::Bool`: (optional) Whether to resume reading from a previous checkpoint. Default is `true`.

# Returns
- The walker object read from the file.

"""
function read_single_walker(filename::String; pbc::String="TTT", resume::Bool=true)
    at = read_single_config(filename; pbc=pbc)
    return convert_system_to_walker(at, resume)    
end

"""
    read_walkers(filename::String; pbc::String="TTT", resume::Bool=true)

Reads walker configurations from a file.

# Arguments
- `filename::String`: The name of the file to read the walker configurations from.
- `pbc::String`: A string specifying the periodic boundary conditions. Default is "TTT".
- `resume::Bool`: A boolean indicating whether to resume reading from a previous checkpoint. Default is `true`.

# Returns
An array of walker objects.

"""
function read_walkers(filename::String; pbc::String="TTT", resume::Bool=true)
    ats = read_configs(filename; pbc=pbc)
    return convert_system_to_walker.(ats, resume)
end

"""
    convert_walker_to_system(at::AtomWalker)

Converts an `AtomWalker` object to an `AbstractSystem` object.

# Arguments
- `at::AtomWalker`: The `AtomWalker` object to be converted.

# Returns
- `AbstractSystem`: The converted `AbstractSystem` object.

"""
function convert_walker_to_system(at::AtomWalker)
    config::FastSystem = at.configuration
    energy = at.energy.val
    iter = at.iter
    num = join(at.list_num_par, ",")
    frozen = join(at.frozen, ",")
    e_frozen = at.energy_frozen_part.val
    return AbstractSystem(config; energy=energy, iter=iter, list_num_par=num, frozen=frozen, energy_frozen_part=e_frozen)
end

"""Convert an AtomicLattice walker to an extXYZ-compatible atomic system."""
function convert_walker_to_system(at::LatticeWalker)
    lattice = at.configuration
    lattice isa AtomicLattice || throw(ArgumentError(
        "only AtomicLattice walkers have an atomic-system representation"))
    sync_ase_lattice!(lattice)
    ase_config = pyconvert(AbstractSystem, lattice.ase_lattice)
    # Strip ASE-only metadata such as `adsorbate_info` (a Python dictionary),
    # which ExtXYZ cannot represent, while preserving atoms, cell and PBC.
    config = FastSystem(
        [Atom(atomic_symbol(atom), position(atom)) for atom in ase_config.particles],
        cell_vectors(ase_config), periodicity(ase_config))
    geometry = string(nameof(typeof(lattice).parameters[2]))
    metadata = (;
        atomic_lattice_schema=1,
        energy=at.energy.val,
        iter=at.iter,
        energy_convention="bare_E_v1",
        freebird_walker="AtomicLattice",
        lattice_geometry=geometry,
        lattice_atom=_encode_atomic_symbols((lattice.lattice_atom,)),
        ase_surface=string(lattice.surface),
        adsorbate_atoms=_encode_atomic_symbols(lattice.adsorbate_atoms),
        supercell_dimensions=join(lattice.supercell_dimensions, ','),
        lattice_constant=_encode_atomic_geometry_scalar(lattice.lattice_constant),
        lattice_periodicity=join(Int.(lattice.periodicity), ','),
        num_nearest_neighbors=lattice.num_nearest_neighbors,
        image_multiplicity=Int(lattice.image_multiplicity),
        type_of_sites=join(lattice.type_of_sites, ','),
        adsorbate_height=_encode_atomic_geometry_scalar(lattice.adsorbate_height),
        site_geometry_fingerprint=_atomic_site_geometry_fingerprint(lattice),
        # The non-numeric prefix prevents ExtXYZ from parsing a long bit mask
        # as a fixed-width integer (and overflowing for realistic lattices).
        component_occupations=join(
            ("bits=" * join(Int.(component)) for component in lattice.components), ';'))
    if lattice.lattice_constant_c !== nothing
        metadata = merge(metadata, (lattice_constant_c=
            _encode_atomic_geometry_scalar(lattice.lattice_constant_c),))
    end
    return AbstractSystem(config; metadata...)
end

"""Identify the walker representation already stored in a trajectory file."""
function _trajectory_storage_kind(filename::String)
    (!isfile(filename) || filesize(filename) == 0) && return nothing
    return open(filename, "r") do io
        first_line = strip(readline(io))
        if tryparse(Int, first_line) !== nothing
            eof(io) && return :unknown
            comment = readline(io)
            occursin(r"(?:^|\s)[\"']?(?:atomic_lattice_schema|freebird_walker)[\"']?\s*=", comment) ||
                return :unknown
            # Parse the frame so quoted values and whitespace around '='
            # follow the same rules as the trajectory reader.
            frames = ExtXYZ.read_frames(filename, 1:1)
            isempty(frames) && return :unknown
            metadata = frames[1]["info"]
            atomic = haskey(metadata, "atomic_lattice_schema") ||
                     get(metadata, "freebird_walker", "") == "AtomicLattice"
            return atomic ?
                   :atomic_lattice_extxyz : :atom_extxyz
        end
        return :unknown
    end
end

function _require_append_storage(filename::String, expected::Symbol)
    existing = _trajectory_storage_kind(filename)
    expected !== :atomic_lattice_extxyz && existing !== :atomic_lattice_extxyz &&
        return nothing
    (existing === nothing || existing === expected) && return nothing
    throw(ArgumentError(
        "cannot append $expected data to '$filename', which contains $existing data"))
end

"""
    append_walker(filename::String, at::AtomWalker)

Append an `AtomWalker` object to a file.
"""
function append_walker(filename::String, at::AtomWalker{C}) where C
    _require_append_storage(filename, :atom_extxyz)
    sys = convert_walker_to_system(at)
    write_frame(filename, ExtXYZ.write_dict(Atoms(sys)), append=true)
end

"""
    append_system(ats1::FlexibleSystem, ats2::FlexibleSystem)

Append two `FlexibleSystem` objects into a single `FastSystem` object.
The first argument is the system to be appended to, and its bounding box and boundary conditions 
will be used for the new system.

# Arguments
- `ats1::FlexibleSystem`: The base system to be appended.
- `ats2::FlexibleSystem`: The system to append.

# Returns
- `new_list`: A new `FastSystem` object containing the appended systems.

"""
function append_system(ats1::FlexibleSystem, ats2::FlexibleSystem)
    new_list = [Atom(atomic_symbol(i),position(i)) for i in ats1.particles]
    append!(new_list,[Atom(atomic_symbol(i),position(i)) for i in ats2.particles])
    return FastSystem(new_list, cell_vectors(ats1), periodicity(ats1))
end

"""
    append_system(ats1::FlexibleSystem, ats2::Vector{FlexibleSystem})
Append a `FlexibleSystem` object to a vector of `FlexibleSystem` objects.
The first argument is the system to be appended to, and its bounding box and boundary conditions
will be used for the new systems.
# Arguments
- `ats1::FlexibleSystem`: The base system to be appended.
- `ats2::Vector{FlexibleSystem}`: A vector of `FlexibleSystem` objects to append.
# Returns
- `configs`: A vector of `FastSystem` objects containing the appended systems.
"""
function append_system(ats1::FlexibleSystem, ats2::Vector{T}) where T
    configs = Vector{FastSystem}(undef, length(ats2))
    Threads.@threads for i in eachindex(ats2)
        configs[i] = append_system(ats1, ats2[i])
    end
    return configs
end

"""
    write_single_walker(filename::String, at::AtomWalker)

Write a single AtomWalker object to a file.

# Arguments
- `filename::String`: The name of the file to write to.
- `at::AtomWalker`: The AtomWalker object to write.

"""
function write_single_walker(filename::String, at::AtomWalker)
    flex = convert_walker_to_system(at)
    save_system(filename::String, flex)
end

function write_single_walker(filename::String, at::LatticeWalker)
    write_walkers(filename, [at])
end

"""
    write_walkers(filename::String, ats::Vector{AtomWalker})

Write a collection of `AtomWalker` objects to a file.

# Arguments
- `filename::String`: The name of the file to write the walkers to.
- `ats::Vector{AtomWalker}`: The collection of `AtomWalker` objects to write.

"""
function write_walkers(filename::String, ats::Vector{AtomWalker{C}}) where C
    flex = convert_walker_to_system.(ats)
    save_trajectory(filename::String, flex)
end

"""
    write_walkers(filename::String, ats::Vector{LatticeWalker})

Write a collection of `LatticeWalker` objects to a file.

# Arguments
- `filename::String`: The name of the file to write the walkers to.
- `ats::Vector{LatticeWalker}`: The collection of `LatticeWalker` objects to write.

"""
function write_walkers(filename::String, ats::Vector{LatticeWalker{C}}; append=false) where C
    # An empty append cannot alter the representation already on disk.
    append && isempty(ats) && isfile(filename) && filesize(filename) > 0 &&
        return nothing
    if !isempty(ats) && all(at -> at.configuration isa AtomicLattice, ats)
        append && _require_append_storage(filename, :atomic_lattice_extxyz)
        systems = convert_walker_to_system.(ats)
        if append
            for system in systems
                write_frame(filename, ExtXYZ.write_dict(Atoms(system)), append=true)
            end
        else
            save_trajectory(filename, systems)
        end
        return nothing
    elseif any(at -> at.configuration isa AtomicLattice, ats)
        throw(ArgumentError("cannot serialize mixed AtomicLattice and non-atomic lattice walkers"))
    end
    append && _require_append_storage(filename, :lattice_csv)
    occupancies = [at.configuration.components for at in ats]
    energies = [at.energy.val for at in ats]
    iter = [at.iter for at in ats]
    df = DataFrame(iter=iter, energy=energies, configuration=occupancies)
    CSV.write(filename, df; append=append)
end

function append_walker(filename::String, at::LatticeWalker{C}) where C
    write_walkers(filename, [at]; append=true)
end

"""
    write_single_walker(filename::String, at::AtomWalker, append::Bool)

Write a single AtomWalker object to a file. If the file already exists, append the walker to the file.

# Arguments
- `filename::String`: The name of the file to write to.
- `at::AtomWalker`: The AtomWalker object to write.
- `append::Bool`: A boolean indicating whether to append the walker to the file if it already exists.

"""
function write_single_walker(filename::String, at::AbstractWalker, append::Bool)
    if isfile(filename) && append
        append_walker(filename, at)
    else
        write_walkers(filename, [at])
    end
end

"""
    extract_free_par(walker::AtomWalker)

Extract free particles from existing walker and create new walker containing only the free particles.

# Arguments
- `walker::AtomWalker`: The AtomWalker object for extraction.

"""

function extract_free_par(walker::AtomWalker)
    free_part = []
    free_indices = []
    components = split_components(walker.configuration, walker.list_num_par)
    for ind in eachindex(components)
        if !walker.frozen[ind]
            push!(free_indices, length(components[ind]))
            for i in 1:length(components[ind])
                pos = components[ind].position[i]
                sym = atomic_symbol(components[ind], i)
                push!(free_part, Atom(sym, pos))
            end
        end
    end
    fast = FastSystem(free_part, cell_vectors(walker.configuration), periodicity(walker.configuration))
    return AtomWalker{length(free_indices)}(fast;list_num_par = free_indices, energy = walker.energy - walker.energy_frozen_part, iter = walker.iter)
end

"""
    generate_random_starting_config(volume_per_particle::Float64, num_particle::Int;
                                    particle_type::Symbol=:H,
                                    periodicity::NTuple{3,Bool}=(false, false, false),
                                    min_separation::Float64=0.0,
                                    max_attempts::Int=1000)

Generate a random starting configuration for a system of particles.

Intended for Metropolis and grand-canonical Metropolis starting states and
reference-run inputs. Nested-sampling live sets must remain i.i.d. draws from the
sampling prior: a minimum-separation draw is not the uniform prior, and a
non-prior initial live set biases the nested-sampling evidence (see the
initialization note in `test/test-SamplingSchemes/test-atomistic-gcns-fixed-n.jl`).

# Arguments
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Int`: The number of particles.
- `particle_type::Symbol=:H`: The type of particle (default is hydrogen).
- `periodicity::NTuple{3,Bool}=(false, false, false)`: Per-axis boundary
  conditions of the generated system. The default keeps the historical
  non-periodic behavior.
- `min_separation::Float64=0.0`: Minimum pairwise separation in Å, enforced by
  sequential insertion with per-particle retries. Distances are minimum-image on
  periodic axes and plain Euclidean otherwise. The default performs no
  separation screening and consumes exactly the historical random stream.
- `max_attempts::Int=1000`: Retry budget per particle; exhaustion throws an
  `ArgumentError` reporting the attempted packing fraction.

# Returns
- `FastSystem`: A FastSystem object representing the generated system.

"""
function generate_random_starting_config(volume_per_particle::Float64, num_particle::Int;
                                         particle_type::Symbol=:H,
                                         periodicity::NTuple{3,Bool}=(false, false, false),
                                         min_separation::Float64=0.0,
                                         max_attempts::Int=1000)
    # generate random starting configuration
    total_volume = volume_per_particle * num_particle
    box_length = total_volume^(1/3)
    box = [[box_length, 0.0, 0.0], [0.0, box_length, 0.0], [0.0, 0.0, box_length]]u"Å"
    boundary_conditions = periodicity
    placed = Vector{Vector{typeof(1.0u"Å")}}()
    for _ in 1:num_particle
        attempts = 0
        while true
            # Literal legacy draw: at min_separation = 0 the first candidate
            # always passes, so the default path consumes exactly the
            # historical stream and returns bit-identical configurations
            candidate = [rand(), rand(), rand()] .* box_length * u"Å"
            if min_separation <= 0.0 ||
               _separation_ok(candidate, placed, box_length, periodicity, min_separation)
                push!(placed, candidate)
                break
            end
            attempts += 1
            if attempts >= max_attempts
                packing = (length(placed) + 1) * (π / 6) * min_separation^3 / total_volume
                throw(ArgumentError(
                    "sequential insertion failed after $max_attempts attempts at particle $(length(placed) + 1) of $num_particle " *
                    "(attempted packing fraction $(round(packing; sigdigits=3)) at the $min_separation Å separation); " *
                    "reduce min_separation or use generate_lattice_starting_config"))
            end
        end
    end
    list_of_atoms = [particle_type => placed[i] for i in 1:num_particle]
    system = atomic_system(list_of_atoms, box, boundary_conditions)
    return FastSystem(system)
end

# Pairwise separation screen for sequential insertion: minimum-image on
# periodic axes, plain Euclidean otherwise
function _separation_ok(candidate, placed, box_length::Float64,
                        periodicity::NTuple{3,Bool}, min_separation::Float64)
    for p in placed
        d2 = 0.0
        for k in 1:3
            δ = ustrip(u"Å", candidate[k] - p[k])
            if periodicity[k]
                δ -= box_length * round(δ / box_length)
            end
            d2 += δ^2
        end
        sqrt(d2) < min_separation && return false
    end
    return true
end

"""
    generate_lattice_starting_config(box_length::Float64, num_particle::Int;
                                     particle_type::Symbol=:H,
                                     periodicity::NTuple{3,Bool}=(true, true, true),
                                     jitter::Float64=0.0)

Generate a starting configuration on a simple-cubic grid: the first
`num_particle` cell-centered sites of the smallest grid containing them, with an
optional uniform jitter of up to `jitter` Å per coordinate. Deterministic at
`jitter = 0`. Covers densities beyond the sequential-insertion regime of
`generate_random_starting_config` (dense-liquid and solid-like starts for
Metropolis annealing); the pairwise separation is at least the grid spacing
minus `2 √3 jitter`. The same nested-sampling scoping applies: this is not a
prior draw and must not seed a nested-sampling live set.

# Arguments
- `box_length::Float64`: The cubic box edge in Å.
- `num_particle::Int`: The number of particles.
- `particle_type::Symbol=:H`: The type of particle.
- `periodicity::NTuple{3,Bool}=(true, true, true)`: Per-axis boundary conditions.
- `jitter::Float64=0.0`: Uniform per-coordinate displacement bound in Å.

# Returns
- `FastSystem`: A FastSystem object representing the generated system.

"""
function generate_lattice_starting_config(box_length::Float64, num_particle::Int;
                                          particle_type::Symbol=:H,
                                          periodicity::NTuple{3,Bool}=(true, true, true),
                                          jitter::Float64=0.0)
    num_particle > 0 || throw(ArgumentError("num_particle must be positive"))
    box = [[box_length, 0.0, 0.0], [0.0, box_length, 0.0], [0.0, 0.0, box_length]]u"Å"
    n_side = ceil(Int, cbrt(num_particle))
    a = box_length / n_side
    list_of_atoms = []
    count = 0
    for ix in 0:(n_side - 1), iy in 0:(n_side - 1), iz in 0:(n_side - 1)
        count >= num_particle && break
        pos = [(ix + 0.5) * a, (iy + 0.5) * a, (iz + 0.5) * a]
        if jitter > 0.0
            pos = pos .+ jitter .* (2 .* [rand(), rand(), rand()] .- 1)
        end
        push!(list_of_atoms, particle_type => pos .* u"Å")
        count += 1
    end
    system = atomic_system(list_of_atoms, box, periodicity)
    return FastSystem(system)
end

"""
    generate_multi_type_random_starting_config(volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])

Generate a random starting configuration for a system of particles with multiple types.

# Arguments
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Vector{Int}`: The number of particles of each type.
- `particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)]`: The types of particles.

# Returns
- `FastSystem`: A FastSystem object representing the generated system.

"""
function generate_multi_type_random_starting_config(volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])
    num_types = length(num_particle)
    total_num_particle = sum(num_particle)
    total_volume = volume_per_particle * total_num_particle
    box_length = total_volume^(1/3)
    box = [[box_length, 0.0, 0.0], [0.0, box_length, 0.0], [0.0, 0.0, box_length]]u"Å"
    boundary_conditions = (false, false, false)
    list_of_atoms = []
    for i in 1:num_types
        for _ in 1:num_particle[i]
            push!(list_of_atoms, particle_types[i] => [rand(), rand(), rand()] .* box_length * u"Å")
        end
    end
    flex = atomic_system(list_of_atoms, box, boundary_conditions)
    return FastSystem(flex)
end

"""
    generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)

Generate initial configurations for a given number of walkers.

# Arguments
- `num_walkers::Int`: The number of walkers.
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Int`: The number of particles.
- `particle_type::Symbol=:H`: The type of particle (default is :H).

# Returns
An array of initial configurations for each walker.

"""
function generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Int; particle_type::Symbol=:H)
    [generate_random_starting_config(volume_per_particle, num_particle; particle_type=particle_type) for _ in 1:num_walkers]
end

"""
    generate_accessible_configs(num_walkers::Int, field::AbstractExternalField, box, num_particle::Int;
                                periodicity::NTuple{3,Bool}=(true, true, true), particle_type::Symbol=:H)

Draw `num_walkers` configurations of `num_particle` atoms with every position uniform in
the accessible region of `field` inside the orthorhombic cell `box` (three cell vectors,
lengths). Each atom is drawn by rejection from the cell: three uniforms per attempt,
scaled by the cell lengths, redrawn until `accessible(field, pos, cell)`. With a field
whose region is the whole cell every attempt succeeds, so each atom consumes exactly
three draws. The configurations are i.i.d. draws from the uniform prior on the
accessible region, the prior that canonical nested sampling under an
`ExternalFieldPotential` requires; pass `V = accessible_volume(field, cell)` to
`gc_thermodynamic_stats_fixed_N` when stitching such runs.

# Returns
- `Vector{FastSystem}`: the configurations.
"""
function generate_accessible_configs(num_walkers::Int, field::AbstractExternalField, box,
                                     num_particle::Int;
                                     periodicity::NTuple{3,Bool}=(true, true, true),
                                     particle_type::Symbol=:H)
    num_particle >= 1 || throw(ArgumentError("num_particle must be at least 1"))
    cell = FastSystem(atomic_system([particle_type => [0.0, 0.0, 0.0]u"Å"], box, periodicity))
    cellv = cell_vectors(cell)
    L = (cellv[1][1], cellv[2][2], cellv[3][3])
    ustrip(u"Å^3", accessible_volume(field, cell)) > 0.0 ||
        throw(ArgumentError("the field's accessible region has zero volume in the cell"))
    configs = FastSystem[]
    for _ in 1:num_walkers
        placed = Vector{typeof(1.0u"Å")}[]
        for _ in 1:num_particle
            pos = [rand() * L[1], rand() * L[2], rand() * L[3]]
            while !accessible(field, pos, cell)
                pos = [rand() * L[1], rand() * L[2], rand() * L[3]]
            end
            push!(placed, pos)
        end
        list_of_atoms = [particle_type => placed[i] for i in 1:num_particle]
        push!(configs, FastSystem(atomic_system(list_of_atoms, box, periodicity)))
    end
    return configs
end


"""
    generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])
Generate initial configurations for a given number of walkers with multiple particle types.
# Arguments
- `num_walkers::Int`: The number of walkers.
- `volume_per_particle::Float64`: The volume per particle.
- `num_particle::Vector{Int}`: A vector containing the number of particles of each type.
- `particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)]`: A vector of symbols representing the types of particles (default is hydrogen and oxygen).
# Returns
An array of initial configurations for each walker, where each configuration contains particles of the specified types.
"""
function generate_initial_configs(num_walkers::Int, volume_per_particle::Float64, num_particle::Vector{Int}; particle_types::Vector{Symbol}=[Symbol(:H), Symbol(:O)])
    [generate_multi_type_random_starting_config(volume_per_particle, num_particle; particle_types=particle_types) for _ in 1:num_walkers]
end

end # module
