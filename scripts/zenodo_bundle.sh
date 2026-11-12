#!/usr/bin/env bash
# Builds the Zenodo package of the paper in one command, run on the login node:
#   bash scripts/zenodo_bundle.sh [out dir]        # default /scratch/sadreyev_lab/pschneider/kimonu_zenodo_2026-11
# On the login node it copies the published figures, Data S1 and S2 (vault figures/*.svg, figures/supp/*.svg, *.csv)
# and the manuscript (for the abstract) to <out>/work/, because the job reads scratch only, and submits itself to SLURM
# (partition sadreyevlab, 2 CPUs, 8 GB, log in <out>/work/logs/). The job
#   1. writes the repository at $KIMONU_REV (default HEAD) with the git-ignored frozen records into <out>/work/tree,
#   2. runs every data-reading figure and number script there under figures/record_inputs.py (<out>/work/listing),
#   3. runs make_zenodo_bundle.py, which writes bundle/, upload/, zenodo_metadata.json, description.html, UPLOAD.md,
#   4. unpacks upload/Kimonu_records.zip into <out>/work/test, checks SHA256SUMS, links the records into the unpacked
#      snapshot, reruns the scripts there and compares the source data, Data S1, S2, figs. S1 to S6 and Figure 3 (C, D
#      redrawn) with the archive (figures/svg_same.py), and tables S2, S3 and every number with step 2. The comparison
#      ends the log.
# The published supplementary figures must be named as the scripts write them (FigS1_stage_timing.svg to
# FigS6_float_errors.svg), or make_zenodo_bundle.py stops. VAULT points at another folder of the same layout.
# A rebuild replaces bundle/, upload/ and work/ of <out>. ZENODO_DOI, ARTICLE_DOI, ARCHIVE_VERSION and RELEASE_DATE,
# when set, fill the placeholders of README.md, CITATION.cff and the metadata (make_zenodo_bundle.py). KZ_EXCLUDE lists
# the nodes the job must not use (default node05 to node08, which have a hardware issue since 2026-11-12; KZ_EXCLUDE=
# for none). Slurm 21.08 has no SBATCH_EXCLUDE.
set -u
K=/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py
VAULT=${VAULT:-/home/sadreyev_lab/pschneider/pipelines/Obsidian/Phil/Kimonu}      # read only
OUT=${1:-${OUT:-/scratch/sadreyev_lab/pschneider/kimonu_zenodo_2026-11}}
W=$OUT/work
export KIMONU_REV=$(git -C $K rev-parse --short=7 "${KIMONU_REV:-HEAD}")

if [ -z "${KZ_BUNDLE_JOB:-}" ]; then          # not SLURM_JOB_ID, which an interactive allocation sets too
  rm -rf "$OUT/bundle" "$OUT/upload" "$W"; mkdir -p "$W/figures/supp" "$W/logs"
  cp -p "$VAULT"/figures/*.svg "$W/figures/" && cp -p "$VAULT"/figures/supp/*.svg "$VAULT"/figures/supp/*.csv "$W/figures/supp/"
  cp -p "$VAULT/Kimonu.md" "$W/manuscript.md"
  cp -p "$0" "$W/zenodo_bundle.sh"
  echo "staged $(ls "$W/figures" "$W/figures/supp" | grep -c '\.') figure files; commit $KIMONU_REV"
  EXCL=${KZ_EXCLUDE-node05,node06,node07,node08}
  exec sbatch -J kz_bundle -p sadreyevlab -c 2 --mem=8G -t 06:00:00 -o "$W/logs/bundle_%j.log" ${EXCL:+--exclude=$EXCL} \
       --export=ALL,KZ_BUNDLE_JOB=1,KIMONU_REV=$KIMONU_REV "$W/zenodo_bundle.sh" "$OUT"
fi

SD=$HOME/miniconda3/envs/straindesign/bin/python          # Python 3.12, the environment of the runs
BIO=$HOME/miniconda3/envs/bioinfo/bin/python              # matplotlib 3.8 for the drawing scripts
export PYTHONPATH=$HOME/sd-kernel-nonzeros PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=$W/mplcfg
FROZEN="sweep/cocoa_v110 sweep/cocoa_eff reference_final reference_final_sampled"     # git-ignored parts of the records
DATA="gather_figdata pair_sets seed_spread dataS1_S2 tableS2_S3 paper_numbers"
DRAW="fig4_figS1_S2_benchmark figS3_S4_uptakes_seeds figS5_S6_float fig3AB_hierarchy fig3CD_coverage"
args() { case $1 in dataS1_S2|fig4_figS1_S2_benchmark|figS3_S4_uptakes_seeds|figS5_S6_float|fig3AB_hierarchy) echo "$2" ;;
         tableS2_S3) echo "$2/TableS2.md $2/TableS3.md" ;; fig3CD_coverage) echo "$W/figures/Figure3.svg $2/Figure3.svg" ;; esac; }
run() {   # <python> <script> <listing dir> <out dir> <record_inputs.py>
  $1 "$5" "$3/$2.txt" figures/$2.py $(args $2 $4) > "$4/$2.stdout" 2> "$4/$2.stderr"; echo "exit $? $2 ($(tail -1 "$4/$2.stderr"))"; }

echo "== 1. tree at $KIMONU_REV $(date +%T)"
mkdir -p $W/tree $W/listing $W/out
git -C $K archive $KIMONU_REV | tar -x -C $W/tree
for d in $FROZEN; do rsync -a $K/scripts/$d/ $W/tree/scripts/$d/; done
echo "== 2. record the inputs $(date +%T)"
cd $W/tree/scripts
for s in $DATA; do run $SD $s $W/listing $W/out $W/tree/scripts/figures/record_inputs.py; done
for s in $DRAW; do run $BIO $s $W/listing $W/out $W/tree/scripts/figures/record_inputs.py; done
echo "== 3. build $(date +%T)"
cd $K/scripts
MANUSCRIPT=$W/manuscript.md $SD make_zenodo_bundle.py $OUT $W/figures $W/listing/*.txt || exit 1
echo "== 4. rebuild from the zip alone $(date +%T)"
T=$W/test; mkdir -p $T/unpacked $T/listing $T/out; cd $T/unpacked
unzip -tq $OUT/upload/Kimonu_records.zip && unzip -q $OUT/upload/Kimonu_records.zip || exit 1
A=$T/unpacked/Kimonu_records; cd $A
sha256sum -c --quiet SHA256SUMS && echo "SHA256SUMS verify, $(wc -l < SHA256SUMS) files"
unzip -q Kimonu.py-$KIMONU_REV.zip -d $T/unpacked; C=$T/unpacked/Kimonu.py-$KIMONU_REV
test -d $C/scripts || { echo "no scripts/ in the snapshot"; exit 1; }
LINKS=$(bash $A/link_records.sh $C) || exit 1; eval "$LINKS"
cd $C/scripts
for s in $DATA; do run $SD $s $T/listing $T/out $W/tree/scripts/figures/record_inputs.py; done
for s in $DRAW; do run $BIO $s $T/listing $T/out $W/tree/scripts/figures/record_inputs.py; done
echo "== compare $(date +%T)"
same() { cmp -s "$1" "$2" && echo "same bytes  $3" || $SD -c "import json,sys; a,b=(json.load(open(x)) for x in sys.argv[1:3]); print(('same JSON   ' if a==b else 'DIFFERENT   ')+sys.argv[3])" "$1" "$2" "$3" 2>/dev/null || echo "DIFFERENT   $3"; }
for f in figdata.json figures/pair_counts.json figures/pair_sets.json figures/seed_spread.json; do same $f $A/figures/source_data/$(basename $f) "$f, rebuilt from the zip, against the archive"; done
for f in DataS1_values.csv DataS2_uptake_acr.csv; do same $T/out/$f $A/figures/supp/$f "$f, rebuilt from the zip, against the published file"; done
for f in TableS2.md TableS3.md; do same $T/out/$f $W/out/$f "$f, rebuilt from the zip, against step 2"; done
for f in supp/FigS1_stage_timing supp/FigS2_reference_agreement out/FigS3_uptakes out/FigS4_seeds out/FigS5_certificate out/FigS6_float_errors out/Figure3; do
  g=$(basename $f); case $g in Figure3) p=$A/figures/$g.svg ;; *) p=$A/figures/supp/$g.svg ;; esac
  $BIO figures/svg_same.py $T/$f.svg $p "$g.svg, drawn from the zip, against the published file"; done
diff <(cut -d'|' -f1-4 $T/out/paper_numbers.stdout) <(cut -d'|' -f1-4 $W/out/paper_numbers.stdout) > $T/paper_numbers.diff \
  && echo "same values paper_numbers.py, $(wc -l < $T/out/paper_numbers.stdout) numbers, rebuilt from the zip, against step 2" \
  || { echo "DIFFERENT   paper_numbers.py"; head -20 $T/paper_numbers.diff; }
echo "files read outside the archive: $(cat $T/listing/*.txt | sort -u | grep '^/' | grep -v "^$A/" | grep -v "^$W/figures/" | wc -l)"
cat $T/listing/*.txt | sort -u | grep '^/' | grep -v "^$A/" | grep -v "^$W/figures/" | head
echo "done $(date +%T)"
