#!/usr/bin/env julia
#
# Export the preprocessed network plus COCOA's complex decomposition (A, Y and the complex
# compositions) so that other implementations can be run on exactly the same input.
#
# The MATLAB reference decides concordance from A and the flux space alone. Handing it COCOA's A
# means the three implementations differ only in their decision procedure, not in how they
# enumerate complexes — which is the comparison the correctness question actually needs.
#
# Usage:
#   julia --project=benchmarks/julia benchmarks/export_crn.jl [profile] [model ...]

import COBREXA
import AbstractFBCModels as A
import SBMLFBCModels
using COCOA
using HiGHS
using SparseArrays

include(joinpath(@__DIR__, "prepare_models.jl"))


function triplets(M)
    rows = rowvals(M)
    vals = nonzeros(M)
    out = String[]
    for j in 1:size(M, 2), k in nzrange(M, j)
        push!(out, "[$(rows[k]), $(j), $(repr(vals[k]))]")
    end
    return out
end


function write_crn(path, model)
    Amat, _ = incidence(model; return_ids = true)
    Y, met_ids = complex_stoichiometry(model; return_ids = true)
    lbs, ubs = A.bounds(model)
    S = A.stoichiometry(model)

    size(Amat, 1) == size(Y, 2) ||
        error("incidence rows $(size(Amat, 1)) != complex columns of Y $(size(Y, 2))")

    open(path, "w") do io
        println(io, "{")
        println(io, "  \"n_complexes\": $(size(Amat, 1)),")
        println(io, "  \"n_reactions\": $(size(Amat, 2)),")
        println(io, "  \"n_metabolites\": $(size(Y, 1)),")
        println(io, "  \"reactions\": [", join(["\"$(r)\"" for r in A.reactions(model)], ", "), "],")
        println(io, "  \"metabolites\": [", join(["\"$(m)\"" for m in met_ids], ", "), "],")
        println(io, "  \"lb\": [", join([repr(x) for x in lbs], ", "), "],")
        println(io, "  \"ub\": [", join([repr(x) for x in ubs], ", "), "],")
        println(io, "  \"A\": [", join(triplets(Amat), ", "), "],")
        println(io, "  \"Y\": [", join(triplets(Y), ", "), "],")
        println(io, "  \"S\": [", join(triplets(S), ", "), "]")
        println(io, "}")
    end
end


function main(args)
    profile = isempty(args) ? "concordance" : args[1]
    models = length(args) > 1 ? args[2:end] : DEFAULT_MODELS
    outdir = joinpath(OUTDIR, profile)
    mkpath(outdir)

    for name in models
        src = joinpath(SHARED, name * ".xml")
        isfile(src) || (@warn "missing model, skipping" name src; continue)
        raw = convert(A.CanonicalModel.Model, COBREXA.load_model(src))
        processed = preprocess(raw, profile)
        write_json(processed, joinpath(outdir, name * ".json"))
        dest = joinpath(outdir, name * ".crn.json")
        write_crn(dest, processed)
        @info "written" dest
    end
end

main(ARGS)
