# Print a Godfrey–Isgur meson spectrum as CSV, with its contribution breakdown.
#
#   julia --project=<env with GIModel> spectrum.jl f1 f2 [nmax] [--L SPD] [--solver ho|fd]
#
# f1 is the quark, f2 the antiquark (u d q s c b). Defaults: nmax = 1, L = SPD,
# solver = ho (the paper's certified oscillator method). Masses in GeV.
# The first line is a `#` comment with the numerics provenance.

using GIModel

function parse_args(args)
    positional = String[]
    L = "SPD"
    solver = "ho"
    i = 1
    while i <= length(args)
        if args[i] == "--L"
            L = uppercase(args[i+1]); i += 2
        elseif args[i] == "--solver"
            solver = lowercase(args[i+1]); i += 2
        else
            push!(positional, args[i]); i += 1
        end
    end
    length(positional) in (2, 3) ||
        error("usage: spectrum.jl f1 f2 [nmax] [--L SPD] [--solver ho|fd]")
    nmax = length(positional) == 3 ? parse(Int, positional[3]) : 1
    solver in ("ho", "fd") || error("--solver must be ho or fd")
    (f1 = Symbol(positional[1]), f2 = Symbol(positional[2]), nmax,
        L_labels = Tuple(string.(collect(L))),
        solver = solver == "ho" ? OscillatorSolver() : FiniteDifferenceSolver())
end

function main(args)
    opts = parse_args(args)
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    meson = Meson(mq, opts.f1, opts.f2)
    spec = compute_spectrum(params, meson;
        levels = spectrum_levels(opts.nmax; L_labels = opts.L_labels),
        solver = opts.solver)
    println("# ", numerics_provenance(opts.solver))
    println("label,n,L,S2p1,J,mass_GeV,central_GeV,contact_GeV,fine_structure_GeV,mixed")
    for s in spec.states
        println(join((s.label, s.n, s.L, s.multiplicity, s.J,
            round(s.mass_GeV; digits = 5), round(s.central_GeV; digits = 5),
            round(s.contact_shift_GeV; digits = 5),
            round(s.fine_structure_shift_GeV; digits = 5),
            !isempty(s.mixings)), ","))
    end
end

main(ARGS)
