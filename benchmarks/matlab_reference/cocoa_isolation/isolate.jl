# Feed COCOA's kinetic stage Kimonu's exact concordance partition, so both kinetic stages see the same input.
import COBREXA; import AbstractFBCModels as A; using COCOA, JSON3, SparseArrays
crn = JSON3.read(read(ARGS[1], String), Dict); part = JSON3.read(read(ARGS[2], String), Dict)
nr, nm = crn["n_reactions"], crn["n_metabolites"]; rxn_ids = String.(crn["reactions"]); met_ids = String.(crn["metabolites"])
S = spzeros(nm, nr); for t in crn["S"]; S[Int(t[1]), Int(t[2])] = Float64(t[3]); end
model = A.CanonicalModel.Model()
for m in met_ids; model.metabolites[m] = A.CanonicalModel.Metabolite(); end
for j in 1:nr
    st = Dict{String,Float64}(); for i in findall(!iszero, S[:, j]); st[met_ids[i]] = S[i, j]; end
    model.reactions[rxn_ids[j]] = A.CanonicalModel.Reaction(lower_bound=Float64(crn["lb"][j]), upper_bound=Float64(crn["ub"][j]), stoichiometry=st)
end
mods = Vector{Set{Symbol}}(); push!(mods, Set(Symbol.(String.(part["balanced"]))))
for m in part["modules"]; push!(mods, Set(Symbol.(String.(m)))); end
println("input: balanced ", length(mods[1]), "  modules ", length(mods)-1)
Amat, cids = COCOA.incidence(model; return_ids=true)
Amat, _, cids = COCOA.augment_with_zero_complex(Amat, COCOA.complex_stoichiometry(model; return_ids=true)[1], cids)
n = length(cids); c2i = Dict(id => i for (i, id) in enumerate(cids))
cadj = COCOA.build_cached_adjacency(Amat, n); gterm = COCOA.compute_global_terminal_complexes(cadj, n)
println("COCOA global_terminal complexes: ", length(gterm), " of ", n)
for eff in (true, false)
    r = COCOA.kinetic_analysis(mods, model; efficient=eff)
    km = r.kinetic_modules; ns = [m for m in km if length(m) > 1]
    cov = reduce(∪, ns; init=Set{Symbol}()); tin = count(c -> haskey(c2i, c) && (c2i[c] in gterm), cov)
    println("efficient=", eff, ": non-singleton ", length(ns), "  largest ", maximum(length, km), "  covered ", length(cov), "  globally-terminal inside ", tin)
end
