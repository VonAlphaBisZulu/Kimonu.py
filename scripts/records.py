"""Where the scripts read and write the records of the paper.

Every result lives in one folder, the records, in the layout of the Zenodo archive (ZENODO_README.md, "Layout"):
networks/, results/kimonu/, results/cocoa_jl/, results/reference/, results/hardware/, figures/source_data/ and
figures/. The unpacked archive (Kimonu_records/) is such a folder. Two environment variables place it:

  KIMONU_RECORDS  the records folder. Every script reads its inputs from it, the networks as well as earlier results.
                  Default: records/ at the top of the repository (ignored by git).
  KIMONU_OUT      where the scripts write, in the same layout. Default: KIMONU_RECORDS. Set it to a new folder to
                  rerun a script without overwriting the archive. A script reads from KIMONU_RECORDS only, so to chain
                  new results into later scripts, point KIMONU_RECORDS at a folder that holds them.

    import records
    records.network("iJR904")                     # <KIMONU_RECORDS>/networks/unblocked/iJR904_fixed.crn.json
    records.read(records.TIMING, "det")           # <KIMONU_RECORDS>/results/kimonu/timing/det
    records.write(records.UPTAKE_MODES, "x.json") # <KIMONU_OUT>/results/kimonu/uptake_modes/x.json, folder made

From a shell (sbatch files): `python records.py read TIMING det` or `python records.py write UPTAKE_MODES` prints the
path (and makes the folder for write).
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECORDS = os.path.abspath(os.environ.get("KIMONU_RECORDS") or os.path.join(REPO, "records"))
OUT = os.path.abspath(os.environ.get("KIMONU_OUT") or RECORDS)

# networks: the input of every tool
UNBLOCKED = "networks/unblocked"                        # <model>_fixed.crn.json, blocked reactions removed
WITH_BLOCKED = "networks/with_blocked_reactions"        # <model>_fixed.crn.json, before the removal
COUPLING_MATRICES = "networks/reference_coupling_matrices"
DEFICIENCY = "networks/deficiency"
# Kimonu
TIMING = "results/kimonu/timing"                        # arms det, att, ref, tan, tanatt, small
UPTAKE_MODES = "results/kimonu/uptake_modes"
HIERARCHY_DETACHED = "results/kimonu/hierarchy_detached"
MEMBER_COUNTS = "results/kimonu/member_counts"
FB_VS_MA = "results/kimonu/flux_balance_vs_mass_action"
UPTAKE_LISTS = "results/kimonu/acr_lists_by_uptake_mode"
GAP_ROOTS = "results/kimonu/gap_roots"
SECOND_ORDER = "results/kimonu/second_order"
INJECTIVITY = "results/kimonu/gap_experiments/injectivity"
DESCARTES = "results/kimonu/gap_experiments/descartes"
AUTONOMOUS_BLOCKS = "results/kimonu/gap_experiments/autonomous_blocks"
TWO_POINT = "results/kimonu/gap_experiments/two_point"
UNIFORM_FLUX = "results/kimonu/uniform_flux"
GENERICITY_BOUND = "results/kimonu/genericity_bound"
LP_ROUNDING = "results/kimonu/lp_point_rounding"
FLOAT_PARTS = "results/kimonu/float_comparison/parts"
LOCAL_ACR = "results/kimonu/symbolic_local_acr"
PREPROCESSING = "results/kimonu/preprocessing_time"
ONE_LP = "results/kimonu/one_lp_vs_fva"
UPSTREAM_ROUTE = "results/kimonu/upstream_set_route_2026-09"   # superseded arm, still read by gather_figdata.py
# frozen records of the other tools
COCOA_FALSE = "results/cocoa_jl/efficient_false"
COCOA_TRUE = "results/cocoa_jl/efficient_true"
COCOA_KINETIC_ILJ478 = "results/cocoa_jl/kinetic_stage_iLJ478"
COCOA_ON_KIMONU = "results/cocoa_jl/kinetic_stage_on_kimonu_partition"
COCOA_NOT_FINISHED = "results/cocoa_jl/runs_not_finished"      # the runs of the networks without a result
REF_EXHAUSTIVE = "results/reference/exhaustive"
REF_SAMPLED = "results/reference/sampled"
REF_KINETIC_ILJ478 = "results/reference/kinetic_stage_iLJ478"
REF_GATE = "results/reference/fork_gate_runs"
REF_TRIAL = "results/reference/trial_runs_2026-08"          # 24-h trial runs before the adaptations of the fork
HARDWARE = "results/hardware"
# figures
SOURCE_DATA = "figures/source_data"                     # figdata.json, pair_counts.json, pair_sets.json, ...
FIGURES = "figures"                                     # the published figures, supp/ for figs. S1 to S6, Data S1, S2
DRAFTS = "drafts"                                       # figure drafts the scripts draw (not part of the archive)
# runs that no paper item reads, kept out of the archive (make_zenodo_bundle.py stops on them)
EXPERIMENTS = "experiments"


def _name(part):
    return globals()[part] if part.isupper() and part in globals() else part


def read(*parts):
    """a path in the records (KIMONU_RECORDS)"""
    if not os.path.isdir(RECORDS):
        raise FileNotFoundError(f"records folder {RECORDS} not found: set KIMONU_RECORDS to the unpacked Zenodo archive "
                                "(Kimonu_records/) or to the folder of the runs")
    return os.path.join(RECORDS, *parts)


def write(*parts, folder=False):
    """a path to write (KIMONU_OUT); makes its folder, or the path itself with folder=True"""
    p = os.path.join(OUT, *parts)
    os.makedirs(p if folder else os.path.dirname(p), exist_ok=True)
    return p


def network(model, blocked=False):
    """the crn.json of a network as every tool received it (blocked=True: before blocked-reaction removal)"""
    return read(WITH_BLOCKED if blocked else UNBLOCKED, f"{model}_fixed.crn.json")


def models():
    """the 33 networks, in the order of kinetic_variants/models33.txt"""
    return [l.strip() for l in open(os.path.join(REPO, "scripts", "kinetic_variants", "models33.txt")) if l.strip()]


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("read", "write"):
        sys.exit(__doc__)
    names = [_name(p) for p in sys.argv[2:]]
    print(read(*names) if sys.argv[1] == "read" else write(*names, folder=True))
