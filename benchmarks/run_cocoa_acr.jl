#!/usr/bin/env julia
#
# Run COCOA's ACR/ACRR on a pre-decomposed network, so its answer can be put beside Kimonu's and
# the reference's published files on the same model.
#
# The exhaustive path (kinetic_efficient=false) is used: the fast pairwise path is known to
# under-report, and the reference's published files come from a pairwise rule of the same kind, so
# using the fast path here would compare two incomplete methods and call it agreement.
#
# Usage:
#   julia --project=benchmarks/julia_revision benchmarks/run_cocoa_acr.jl <crn.json> <out.json> [seed]

import COBREXA
import AbstractFBCModels as A
using COCOA
using HiGHS
using Gurobi
if get(ENV, "COCOA_SOLVER", "HiGHS") == "CPLEX"
    using CPLEX
end

const SOLVER = get(ENV, "COCOA_SOLVER", "HiGHS")

# COCOA batches candidate pairs across Distributed workers; passing an empty vector disables that
# entirely, which is how every earlier timing here was produced.
using Distributed
const SOLVER_EARLY = get(ENV, "COCOA_SOLVER", "HiGHS")
const NPROCS = parse(Int, get(ENV, "COCOA_PROCS", "0"))
if NPROCS > 0
    addprocs(NPROCS; exeflags = "--project=$(Base.active_project())")
    @everywhere using COCOA
    @everywhere using HiGHS
    if SOLVER_EARLY == "Gurobi"
        @everywhere using Gurobi
    elseif SOLVER_EARLY == "CPLEX"
        @everywhere using CPLEX
    end
end
const COCOA_WORKERS = NPROCS > 0 ? workers() : Int[]
function optimizer()
    SOLVER == "Gurobi" && return Gurobi.Optimizer
    SOLVER == "CPLEX" && return Base.invokelatest(getfield, Main, :CPLEX).Optimizer
    return HiGHS.Optimizer
end
using SparseArrays
using JSON3

crn_path, out_path = ARGS[1], ARGS[2]
seed = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 1234

crn = JSON3.read(read(crn_path, String), Dict)
nr, nm = crn["n_reactions"], crn["n_metabolites"]
rxn_ids = String.(crn["reactions"])
met_ids = String.(crn["metabolites"])

S = spzeros(nm, nr)
for t in crn["S"]
    S[Int(t[1]), Int(t[2])] = Float64(t[3])
end

model = A.CanonicalModel.Model()
for m in met_ids
    model.metabolites[m] = A.CanonicalModel.Metabolite()
end
for j in 1:nr
    stoich = Dict{String,Float64}()
    for i in findall(!iszero, S[:, j])
        stoich[met_ids[i]] = S[i, j]
    end
    model.reactions[rxn_ids[j]] = A.CanonicalModel.Reaction(
        lower_bound = Float64(crn["lb"][j]), upper_bound = Float64(crn["ub"][j]),
        stoichiometry = stoich)
end

@info "model rebuilt" reactions = length(A.reactions(model)) metabolites = length(A.metabolites(model))

# activity_concordance_analysis does not forward `seed` to its kinetic stage, but that stage is
# deterministic on this path: kinetic_analysis calls identify_acr_acrr, which uses no RNG. The
# seeded identify_acr_acrr_dce is exported and never reached. So the seed varies the concordance
# sampler only, which is the stage under comparison.
t_conc = @elapsed res = activity_concordance_analysis(
    model; optimizer = optimizer(), workers = COCOA_WORKERS, seed = UInt(seed),
    kinetic_analysis = true, kinetic_efficient = false)

acr = collect(res.acr.metabolite_id)
acrr = [[a, b] for (a, b) in zip(res.acrr.metabolite_1, res.acrr.metabolite_2)]
kinetic = collect(res.complexes.kinetic_module)
giant = isempty(kinetic) ? 0 : maximum([count(==(k), kinetic) for k in unique(kinetic) if k > 0];
                                       init = 0)

open(out_path, "w") do io
    JSON3.write(io, Dict("tool" => "COCOA.jl", "seed" => seed, "seconds" => t_conc,
                         # per-complex module labels, not just the giant's size: without them a
                         # disagreement cannot be attributed to the partition rather than to the
                         # span test that reads it.
                         "kinetic_module" => kinetic,
                         "complex_ids" => String.(res.complexes.complex_id),
                         "acr" => acr, "acrr" => acrr, "giant_module" => giant,
                         "n_complexes" => length(kinetic)))
end
@info "written" out_path n_acr = length(acr) n_acrr = length(acrr) giant seconds = round(t_conc, digits=1)
