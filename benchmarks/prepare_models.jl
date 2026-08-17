#!/usr/bin/env julia
#
# Stage 1 of the benchmark: run COCOA.jl's own preprocessing pipeline and dump the resulting
# network to JSON, so that every tool in the comparison consumes the *same* preprocessed model.
#
# Preprocessing parity is not a detail. Complex counts (and therefore module and concordant-pair
# counts) depend on which of normalize_bounds / remove_blocked_reactions / remove_orphans /
# split_into_elementary / split_into_irreversible ran. Comparing tools that preprocessed
# differently compares preprocessing, not concordance.
#
# Usage:
#   julia --project=benchmarks/julia benchmarks/prepare_models.jl [profile] [model ...]
#
# Profiles:
#   concordance  normalize_bounds, remove_blocked_reactions, remove_orphans, split_into_irreversible
#   kinetic      the above plus split_into_elementary (ordered mechanism), as COCOA's README
#                recommends for kinetic-module analysis
#
# Output: benchmarks/processed/<profile>/<model>.json

import COBREXA
import AbstractFBCModels as A
import SBMLFBCModels
using COCOA
using HiGHS
using SparseArrays

const HERE = @__DIR__
const SHARED = get(ENV, "KIMONU_SHARED_MODELS", joinpath(HERE, "..", "..", "shared_models"))
const OUTDIR = joinpath(HERE, "processed")

const DEFAULT_MODELS = ["EnvZ_OmpR", "e_coli_core", "e_coli_core_anaerobic", "iAB_RBC_283"]


"""Write a canonical model as JSON: metabolite ids, and per reaction its stoichiometry and bounds.

Deliberately minimal and explicit rather than a standard exchange format — the point is that the
Python side rebuilds exactly this network, with no format-conversion ambiguity in between.
"""
function write_json(model, path)
    mets = A.metabolites(model)
    rxns = A.reactions(model)
    S = A.stoichiometry(model)
    lbs, ubs = A.bounds(model)

    open(path, "w") do io
        println(io, "{")
        println(io, "  \"metabolites\": [", join(["\"$(m)\"" for m in mets], ", "), "],")
        println(io, "  \"reactions\": [")
        rows = rowvals(S)
        vals = nonzeros(S)
        for (j, r) in enumerate(rxns)
            terms = ["\"$(mets[rows[k]])\": $(repr(vals[k]))" for k in nzrange(S, j)]
            comma = j == length(rxns) ? "" : ","
            println(io, "    {\"id\": \"$(r)\", \"lb\": $(repr(lbs[j])), \"ub\": $(repr(ubs[j])), ",
                        "\"stoichiometry\": {", join(terms, ", "), "}}", comma)
        end
        println(io, "  ]")
        println(io, "}")
    end
end


function preprocess(model, profile)
    m = model |>
        normalize_bounds |>
        m -> remove_blocked_reactions(m; optimizer = HiGHS.Optimizer) |>
        remove_orphans
    if profile == "kinetic"
        m = split_into_elementary(m)
    end
    return split_into_irreversible(m)
end


function main(args)
    profile = isempty(args) ? "concordance" : args[1]
    profile in ("concordance", "kinetic") || error("unknown profile: $profile")
    models = length(args) > 1 ? args[2:end] : DEFAULT_MODELS

    outdir = joinpath(OUTDIR, profile)
    mkpath(outdir)

    for name in models
        src = joinpath(SHARED, name * ".xml")
        isfile(src) || (@warn "missing model, skipping" name src; continue)

        @info "preprocessing" name profile
        raw = convert(A.CanonicalModel.Model, COBREXA.load_model(src))
        t = @elapsed processed = preprocess(raw, profile)

        dest = joinpath(outdir, name * ".json")
        write_json(processed, dest)
        @info "written" dest reactions = length(A.reactions(processed)) metabolites =
            length(A.metabolites(processed)) seconds = round(t, digits = 2)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
