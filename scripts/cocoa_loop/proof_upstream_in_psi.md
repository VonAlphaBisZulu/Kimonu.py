# The iterated ψ route contains the published route (Langary et al. 2025, COCOA.jl)

Setting. Network without blocked reactions; A (complexes × reactions, ∅ included as a complex with zero
composition), Y, N = Y·A, S the substrate incidence (S_{c,j} = 1 iff c is the substrate of reaction j).
Kimonu draws r* > 0 in ker N, sets W = N·diag(r*)·Sᵀ and reads coupled sets as proportional rows of a
basis of ker W (u = ψ/ψ*). L := A·diag(r*)·Sᵀ is a Laplacian of the complex graph with positive weights
r* (columns sum to zero, off-diagonals ≥ 0), and W = Y·L. Concordance is Kimonu's exact one: rows of
M = A·K, K = ker N.

Two observations used throughout.
* (O1) u = 1 lies in ker W, and in ker W_E for every set E of equalities u_c = u_d: W·1 = N r* = 0.
* (O2) For every u ∈ ker W_E, r(u) := diag(r*)·Sᵀ·u lies in ker N (N r(u) = W u = 0), so the activity
  a(u) = A r(u) = L u obeys every linear relation that holds on A·ker N: balanced complexes have a_b(u) = 0,
  and a concordance module m has a_m(u) = α_m(u)·μ_m with a fixed vector μ_m and a linear functional α_m.
  No positivity of u is needed, because Kimonu's concordance is a statement about the whole of ker N.

## (a) Upstream sets are proportional rows of ker W — PROVEN

Claim. Let m be a concordance module (exact), C̄m = C_b ∪ m its extended module, and C0 ⊆ C̄m any upstream
set (autonomous and feeding; Langary et al. Supplement S.2.1), in particular the output of the Upstream
Algorithm (Theorem S2-7/S2-9). Then u_c = u_d for all c, d ∈ C0 and all u ∈ ker W. Hence C0 lies in one
class of proportional rows of any basis of ker W.

Proof. Order complexes C0 first. Autonomy ⇒ L is block lower triangular (Lemma S2-1; it uses only the
sparsity pattern, and L has the pattern of A_k with k = r*). Feeding ⇒ L00 is WCDD, hence nonsingular
(Lemma S2-2, Corollary S2-3; again pattern plus positive weights). So u0 = L00⁻¹ a0(u) for every
u ∈ ker W. By (O2), a0(u) = α(u)·(0; μ) on C0 ⊆ C_b ∪ m, so u0 = α(u)·v0 with v0 := L00⁻¹(0; μ) fixed.
By (O1), u = 1 gives 1_{C0} = α(1)·v0, so α(1) ≠ 0 and v0 = 1/α(1)·1: every entry equal. Hence
u_c = α(u)/α(1) for all c ∈ C0. ∎

This is Langary's Proposition S2-6 with ψ replaced by u ∈ ker W. Their proof uses only the linear relation
a = A_k ψ and the concordance structure of a; both hold on all of ker W (O2), so ker W loses nothing.

Uptake modes. COCOA augments the network with ∅ but never puts ∅ into a concordance module, so its upstream
sets never contain ∅ and never contain an uptake product (autonomy). With detached uptakes Kimonu replaces
∅'s column in W by one source per uptake reaction; for a set C0 avoiding ∅ and uptake products, autonomy,
feeding and the block L00 are unchanged, so (a) holds in both modes.

## (b) Every step of the S.4.1 loop is matched at Kimonu's fixed point — PROVEN

Let E∞ be Kimonu's final equality set, V := ker W_{E∞}, and K∞ its final modules. At the fixed point
(F1) the modules of size > 1 are exactly the classes of proportional rows of V (the last round merged
nothing), so u_c = u_d on V within a module (rows proportional and equal at u = 1); (F2) the closure is
stable: no Proposition S3-4 merge and no new ACR metabolite with YΔ(K∞) and Kimonu's known ACR set.

Induction over COCOA's iterations t, with invariants
* (I1) every concordance module of COCOA's current partition P_t is V-concordant (a_m(u) = α_m(u) μ_m on V),
  and balanced complexes have zero activity on V;
* (I2) every current COCOA kinetic module lies inside one module of K∞;
* (I3) COCOA's known ACR set lies inside Kimonu's.

Base. P_0 is Kimonu's exact concordance: V-concordant by (O2).
Upstream step. (a) with V in place of ker W (it needs only (O1), (O2) and I1): every upstream set has
u_c = u_d on V. Upstream sets consist of reactant complexes (Lemma S2-4), so each lies in one V-class of
substrate complexes, i.e. in one module of K∞ (F1).
Overlap merges (Lemma S3-1). Transitive: ✓.
Proposition S3-4 merges (with Remark S3-6). COCOA's YΔ is built from its sets, each inside a module of
K∞, and its ACR augmentation is ⊆ Kimonu's (I3); so im(YΔ_COCOA) ⊆ im(YΔ_Kimonu). A merge licensed for
COCOA is licensed for Kimonu, and by (F2) the two complexes are already in one Kimonu module. I2 ✓.
ACR identification. e_s ∈ im(YΔ_COCOA + known) ⊆ im(YΔ_Kimonu + known): by (F2) Kimonu has it. I3 ✓.
Proposition S4-1 (concordance merge). If upstream sets C1 ⊆ C̄m1, C2 ⊆ C̄m2 end up in one COCOA module,
they lie in one V-class, so u_{C1} = u_{C2} on V. From (a), u_{C1} = α_{m1}(u)/α_{m1}(1) and
u_{C2} = α_{m2}(u)/α_{m2}(1), so α_{m2} = κ·α_{m1} on V: m1 ∪ m2 is V-concordant. I1 ✓.
Balanced linkage classes (add_singleton_balanced, Proposition S2-11). A linkage class whose complexes are
all balanced carries a positive flow (r* > 0) conserved at every node, so it is strongly connected; its
Laplacian block has a one-dimensional kernel, and a_L(u) = L_LL u_L = 0 on V with u = 1 in V gives
u_c = u_d on V for the whole class. ✓
Theorem S4-6 applied to a single unbalanced module. Then all complexes form one extended module whose
upstream set is the set of nonterminal complexes (Theorem S2-7); (a) applies. ✓

Deficiency shortcut (Lemma S4-5, S4-7, Proposition S4-8, Theorem S4-6). COCOA declares δ_k = 1 when
δ − (number of concordance merges) = 1 and the network is not weakly reversible, then merges all unbalanced
concordance modules and couples all nonterminal complexes. The bounds are sound, so the shortcut fires only
when the mass-action deficiency δ_k = dim(ker Y ∩ im A_k) is one. im L = im A_k, since
diag(r*)·Sᵀ = diag(k)·Sᵀ·diag(ψ*). Kimonu's first ker W already holds both conclusions in that case.
ker W = {u : L u ∈ ker Y ∩ im L}, and L·1 = A r* lies in ker Y (N r* = 0) and is nonzero (r* is not
complex-balanced in a network that is not weakly reversible). So ker Y ∩ im L = span(L·1) and
ker W = ker L + span(1). Vectors in ker L vanish on nonterminal complexes (L is a Laplacian with positive
weights; its kernel lives on the terminal strong linkage classes). Hence every u ∈ ker W, and every u in
the smaller ker W_E, has u_c = λ(u) on all nonterminal complexes: they form one class of proportional rows.
Also a(u) = L u = λ(u)·L·1 on ker W, so all unbalanced complexes are V-concordant and the merged module
satisfies I1. ✓ Checked on 300 random deficiency-one networks (202 with two or more nonterminal
complexes): ker W couples all nonterminal complexes in every one (deficiency_one_check.py).

Conclusion. For every network with an exact concordance partition and generic r*, every COCOA kinetic
module (kinetic stage, efficient=false, fed with the same concordance) lies inside one module of Kimonu's
iterated route, and COCOA's ACR set lies inside Kimonu's.

## What the proof does not say

* Strictness: the containment can be equality. Whether Kimonu finds more is an empirical question
  (comparison table), and where it does, it is because ker W relates complexes that are not in one
  upstream set: complexes outside every extended module's upstream set (entry complexes fed from other
  modules), and linear relations across several sets that the proportional-row readout combined with the
  feedback picks up but per-module Laplacian inversion cannot.
* Non-reactant complexes: Kimonu keeps every complex that is no reaction's substrate as a singleton.
  COCOA's final Proposition S3-4 pass can attach such a complex to a module (Langary et al., Fig. 1:
  "all nonterminal complexes along with one terminal complex"). That attachment is valid (ψ ratio fixed by
  the log-space argument alone), so on non-reactant complexes COCOA can be coarser; the containment above
  is about reactant complexes.
