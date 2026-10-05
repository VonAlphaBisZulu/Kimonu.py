#!/usr/bin/env julia
#
# COCOA v1.1.0's own kinetic stage (efficient=false: the full Section S.4.1 loop of Langary et al. --
# Upstream Algorithm, Proposition S3-4 merges, Remark S3-6 ACR augmentation, Proposition S4-1
# concordance merges, Theorem S4-6) run on Kimonu's exact concordance modules instead of COCOA's
# sampled ones. With the concordance input identical, any remaining difference to Kimonu's iterated
# route comes from the kinetic loop itself.
#
# Usage: julia --project=benchmarks/julia_v110 cocoa_side.jl <crn.json> <kimonu_side.json> <out.json>

import AbstractFBCModels as A
using COCOA
using SparseArrays, JSON3, Printf

crn_path, kim_path, out_path = ARGS[1], ARGS[2], ARGS[3]
crn = JSON3.read(read(crn_path, String), Dict)
kim = JSON3.read(read(kim_path, String), Dict)
nr, nm = crn["n_reactions"], crn["n_metabolites"]
rxn_ids = String.(crn["reactions"]); met_ids = String.(crn["metabolites"])

# the same model rebuild as benchmarks/run_cocoa_acr.jl
S = spzeros(nm, nr)
for t in crn["S"]
    S[Int(t[1]), Int(t[2])] = Float64(t[3])
end
model = A.CanonicalModel.Model()
for m in met_ids
    model.metabolites[m] = A.CanonicalModel.Metabolite()
end
for j in 1:nr
    st = Dict{String,Float64}(met_ids[i] => S[i, j] for i in findall(!iszero, S[:, j]))
    model.reactions[rxn_ids[j]] = A.CanonicalModel.Reaction(
        lower_bound = Float64(crn["lb"][j]), upper_bound = Float64(crn["ub"][j]), stoichiometry = st)
end

# COCOA's complex symbols, keyed by composition exactly as kimonu_side.py writes it
Y, ymets, cids = complex_stoichiometry(model; return_ids = true)
Yc = SparseMatrixCSC(Y)
key_of = Dict{Symbol,String}()
sym_of = Dict{String,Symbol}()
for (j, c) in enumerate(cids)
    rows = Yc.rowval[nzrange(Yc, j)]; vals = Yc.nzval[nzrange(Yc, j)]
    parts = sort([(String(ymets[i]), v) for (i, v) in zip(rows, vals)])
    k = join([@sprintf("%s:%.9g", m, v) for (m, v) in parts], "|")
    key_of[c] = k; sym_of[k] = c
end
tosym(ks) = Set{Symbol}(sym_of[k] for k in ks)
missing_keys = [k for k in vcat(String.(kim["balanced"]), String.(reduce(vcat, kim["concordance_modules"]; init=Any[]))) if !haskey(sym_of, k)]
isempty(missing_keys) || error("complexes without a COCOA match: $(length(missing_keys)), e.g. $(missing_keys[1])")

conc = Set{Symbol}[tosym(String.(kim["balanced"]))]
for g in kim["concordance_modules"]
    push!(conc, tosym(String.(g)))
end
@info "input" balanced = length(conc[1]) modules = length(conc) - 1 complexes = length(cids)

t = @elapsed res = kinetic_analysis(conc, model; efficient = false)
mods = [sort([key_of[c] for c in m]) for m in res.kinetic_modules if length(m) > 1]
open(out_path, "w") do io
    JSON3.write(io, Dict("tool" => "COCOA.jl v1.1.0 kinetic stage (efficient=false) on Kimonu concordance",
                         "seconds" => t, "modules" => mods,
                         "acr" => sort(String.(res.acr_metabolites)),
                         "acrr" => sort([sort([String(a), String(b)]) for (a, b) in res.acrr_pairs])))
end
@info "written" out_path seconds = round(t, digits = 1) n_modules = length(mods) n_acr = length(res.acr_metabolites) n_acrr = length(res.acrr_pairs)
