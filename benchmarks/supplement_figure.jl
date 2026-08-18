#!/usr/bin/env julia
#
# Reproduce the three panels of the COCOA supplement (giant kinetic module size, ACR free
# metabolites, ACRR free metabolite pairs) with Kimonu added as a third series.
#
# Kimonu has no kinetic-module or ACR layer of its own yet, so its series is produced by feeding
# Kimonu's exact concordance partition into COCOA's kinetic/ACR stage. That isolates exactly what
# differs between the implementations — the concordance partition — and leaves the downstream
# machinery identical for every series. It is not a pure Kimonu pipeline and must not be labelled
# as one.
#
# COCOA's candidate prefilter is stochastic, so its series is run over several seeds and reported
# with a spread.
#
# Usage:
#   julia --project=benchmarks/julia benchmarks/supplement_figure.jl <crn.json> <kimonu.json> <out.json> [seeds...]

import COBREXA
import AbstractFBCModels as A
using COCOA
using SparseArrays
using JSON3

using HiGHS
using Gurobi

const SOLVER = get(ENV, "COCOA_SOLVER", "HiGHS")

optimizer() = SOLVER == "Gurobi" ? Gurobi.Optimizer : HiGHS.Optimizer

"""Rebuild a canonical model from the shared exported network."""
function model_from_crn(crn)
    nr = crn["n_reactions"]
    nm = crn["n_metabolites"]
    rxn_ids = String.(crn["reactions"])
    met_ids = String.(crn["metabolites"])

    S = spzeros(nm, nr)
    for t in crn["S"]
        S[Int(t[1]), Int(t[2])] = Float64(t[3])
    end

    model = A.CanonicalModel.Model()
    for (i, m) in enumerate(met_ids)
        model.metabolites[m] = A.CanonicalModel.Metabolite()
    end
    for j in 1:nr
        stoich = Dict{String,Float64}()
        for i in findall(!iszero, S[:, j])
            stoich[met_ids[i]] = S[i, j]
        end
        model.reactions[rxn_ids[j]] = A.CanonicalModel.Reaction(
            lower_bound = Float64(crn["lb"][j]),
            upper_bound = Float64(crn["ub"][j]),
            stoichiometry = stoich,
        )
    end
    return model
end

"""Concordance modules as COCOA expects them: sets of complex ids, in the model's complex order."""
function modules_from_labels(labels, complex_ids)
    groups = Dict{Int,Set{Symbol}}()
    for (i, l) in enumerate(labels)
        l <= 0 && continue
        push!(get!(groups, l, Set{Symbol}()), complex_ids[i])
    end
    return collect(values(groups))
end

function panel_values(kin)
    giant = isempty(kin.kinetic_modules) ? 0 : maximum(length.(kin.kinetic_modules))
    return (giant_module = giant,
            n_kinetic_modules = count(m -> length(m) > 1, kin.kinetic_modules),
            n_acr = length(kin.acr),
            n_acrr = length(kin.acrr))
end

function main(args)
    crn_path, kimonu_path, out_path = args[1], args[2], args[3]
    seeds = length(args) > 3 ? parse.(Int, args[4:end]) : [1234, 42, 7, 999, 2026, 31337]

    crn = JSON3.read(read(crn_path, String), Dict)
    model = model_from_crn(crn)
    _, complex_ids = incidence(model; return_ids = true)
    @info "model rebuilt" reactions = length(A.reactions(model)) complexes = length(complex_ids)

    opt = optimizer()
    results = Dict{String,Any}("model" => basename(crn_path), "solver" => SOLVER, "seeds" => seeds)

    # --- COCOA, one run per seed -----------------------------------------------------------
    # Wall time is recorded per stage: the concordance sweep is what differs between the
    # implementations, while the kinetic/ACR stage is shared, so a combined figure would blur the
    # comparison the runtime panel is meant to make.
    cocoa_runs = []
    for s in seeds
        @info "COCOA" seed = s
        t_conc = @elapsed res = activity_concordance_analysis(
            model; optimizer = opt, workers = Int[], seed = UInt(s),
            kinetic_analysis = false)
        t_kin = @elapsed kin = kinetic_analysis(res, model; efficient = false)
        push!(cocoa_runs, merge(Dict(pairs(panel_values(kin))),
                                Dict("seconds_concordance" => t_conc,
                                     "seconds_kinetic" => t_kin,
                                     "seconds" => t_conc + t_kin,
                                     "seed" => s)))
    end
    results["cocoa"] = cocoa_runs

    # --- Kimonu's partition through the same kinetic stage ---------------------------------
    kim = JSON3.read(read(kimonu_path, String), Dict)
    labels = [Int(c["module"]) for c in kim["complexes"]]
    length(labels) == length(complex_ids) ||
        error("Kimonu reports $(length(labels)) complexes, model has $(length(complex_ids))")
    t = @elapsed kin = kinetic_analysis(modules_from_labels(labels, complex_ids), model;
                                        efficient = false)
    results["kimonu"] = merge(Dict(pairs(panel_values(kin))), Dict("seconds" => t))

    open(out_path, "w") do io
        JSON3.write(io, results)
    end
    @info "written" out_path
end

main(ARGS)
