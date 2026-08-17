#!/usr/bin/env julia
#
# Stage 2a of the benchmark: run COCOA.jl's LP-sweep concordance analysis and record both its
# timing and its partition.
#
# The preprocessed model is written to JSON in the same run, so the Python side analyses the exact
# network COCOA analysed. Complexes are exported by *composition* (a metabolite -> coefficient map)
# rather than by id, so the two implementations can be matched without relying on either one's
# naming scheme.
#
# Usage:
#   julia --project=benchmarks/julia benchmarks/run_cocoa.jl [--profile P] [--scenario i|ii]
#                                                            [--seed N] [--repeats N] [--model M]...
#
# Scenarios follow the COCOA paper: (i) is the unbounded steady-state cone, (ii) additionally
# constrains the objective to within a relative tolerance of its optimum.
#
# Output: benchmarks/results/cocoa_<profile>_<scenario>_<model>.json
#         benchmarks/processed/<profile>/<model>.json

import COBREXA
import AbstractFBCModels as A
import SBMLFBCModels
using COCOA
using HiGHS
using SparseArrays

include(joinpath(@__DIR__, "prepare_models.jl"))

const RESULTS = joinpath(@__DIR__, "results")


function parse_args(args)
    opts = Dict{String,Any}("profile" => "concordance", "scenario" => "i",
                            "seed" => 1234, "repeats" => 2, "models" => String[])
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--profile";       opts["profile"] = args[i+1];            i += 2
        elseif a == "--scenario";  opts["scenario"] = args[i+1];           i += 2
        elseif a == "--seed";      opts["seed"] = parse(Int, args[i+1]);   i += 2
        elseif a == "--repeats";   opts["repeats"] = parse(Int, args[i+1]); i += 2
        elseif a == "--model";     push!(opts["models"], args[i+1]);       i += 2
        else error("unknown argument: $a")
        end
    end
    isempty(opts["models"]) && (opts["models"] = DEFAULT_MODELS)
    opts["scenario"] in ("i", "ii") || error("scenario must be i or ii")
    return opts
end


"""Complex compositions, in COCOA's own complex order, as a vector of (met_id, coeff) pairs.

`complex_stoichiometry` returns Y (metabolites x complexes) and the *metabolite* ids — complexes
themselves are identified only by column position, which is exactly what we want: the composition
is the identity used to match against the Python side.
"""
function complex_compositions(model)
    Y, met_ids = complex_stoichiometry(model; return_ids = true)
    rows = rowvals(Y)
    vals = nonzeros(Y)
    compositions = [[(string(met_ids[rows[k]]), vals[k]) for k in nzrange(Y, j)]
                    for j in 1:size(Y, 2)]
    return compositions, ["complex_$(j)" for j in 1:size(Y, 2)]
end


function write_results(path, opts, model_name, compositions, complex_ids, modules, seconds, extra)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"tool\": \"COCOA.jl\",")
        println(io, "  \"model\": \"$(model_name)\",")
        println(io, "  \"profile\": \"$(opts["profile"])\",")
        println(io, "  \"scenario\": \"$(opts["scenario"])\",")
        println(io, "  \"seed\": $(opts["seed"]),")
        println(io, "  \"seconds\": $(repr(seconds)),")
        for (k, v) in extra
            println(io, "  \"$(k)\": $(repr(v)),")
        end
        println(io, "  \"complexes\": [")
        for j in 1:length(complex_ids)
            terms = ["\"$(m)\": $(repr(c))" for (m, c) in compositions[j]]
            comma = j == length(complex_ids) ? "" : ","
            println(io, "    {\"id\": \"$(complex_ids[j])\", \"module\": $(modules[j]), ",
                        "\"composition\": {", join(terms, ", "), "}}", comma)
        end
        println(io, "  ]")
        println(io, "}")
    end
end


function main(args)
    opts = parse_args(args)
    mkpath(RESULTS)
    outdir = joinpath(OUTDIR, opts["profile"])
    mkpath(outdir)

    for name in opts["models"]
        src = joinpath(SHARED, name * ".xml")
        isfile(src) || (@warn "missing model, skipping" name src; continue)

        raw = convert(A.CanonicalModel.Model, COBREXA.load_model(src))
        processed = preprocess(raw, opts["profile"])
        write_json(processed, joinpath(outdir, name * ".json"))

        compositions, complex_ids = complex_compositions(processed)

        # Scenario (ii) restricts the flux space to near-optimal states; (i) is the full cone.
        bound = opts["scenario"] == "ii" ? COBREXA.relative_tolerance_bound(0.999) : nothing

        @info "running COCOA" name profile = opts["profile"] scenario = opts["scenario"] repeats =
            opts["repeats"]
        # The first call pays Julia's JIT compilation; reporting it as the runtime would flatter the
        # Python side. Repeat and report the fastest run, keeping the cold time for the record.
        times = Float64[]
        local result
        for _ in 1:opts["repeats"]
            push!(times, @elapsed result = activity_concordance_analysis(
                processed;
                optimizer = HiGHS.Optimizer,
                objective_bound = bound,
                workers = Int[],                # in-process: the timing is single-threaded
                seed = UInt(opts["seed"]),      # the default is rand(UInt) — pin it or runs differ
                kinetic_analysis = false,
            ))
        end
        seconds = minimum(times)

        modules = collect(result.complexes.concordance_module)
        length(modules) == length(complex_ids) ||
            error("complex count mismatch: analysis $(length(modules)) vs export $(length(complex_ids))")

        dest = joinpath(RESULTS, "cocoa_$(opts["profile"])_$(opts["scenario"])_$(name).json")
        write_results(dest, opts, name, compositions, complex_ids, modules, seconds,
                      ["n_reactions" => length(A.reactions(processed)),
                       "n_metabolites" => length(A.metabolites(processed)),
                       "seconds_cold" => first(times), "repeats" => opts["repeats"]])
        @info "written" dest seconds = round(seconds, digits = 2) complexes = length(complex_ids)
    end
end

main(ARGS)
