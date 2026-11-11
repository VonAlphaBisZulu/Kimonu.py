"""Run write_unblocked_crn.py (fva_time.sbatch) with cobra's process count at 1, so straindesign.fva stays in one process."""
import runpy, sys
import cobra
cobra.Configuration().processes = 1
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name="__main__")
