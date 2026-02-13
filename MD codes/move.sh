#!/usr/bin/env bash
# Copy restarts & slurm_run.sh, then submit all leaf jobs immediately (no wait).
# Run from a *composition* directory that has 1.glass/{restart.*,slurm_run.sh}

SRC_DIR="1.glass"
REQ_FILES=( "restart.glass" "restart.anneal_500K" "restart.anneal_700K" "slurm_run.sh" )
#TARGET_ROOTS=( "2.indent" "3.srjt" "4.creep" )
TARGET_ROOTS=( "3.srjt" "4.creep" )

# --- sanity checks ---
if [[ ! -d "$SRC_DIR" ]]; then
  echo "ERROR: '$SRC_DIR' directory not found. Run this from a composition directory." >&2
  exit 1
fi
for f in "${REQ_FILES[@]}"; do
  if [[ ! -f "$SRC_DIR/$f" ]]; then
    echo "ERROR: missing '$SRC_DIR/$f'." >&2
    exit 1
  fi
done
if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: 'sbatch' not found in PATH." >&2
  exit 1
fi

# --- mapping: state dir -> restart filename ---
restart_for_state() {
  case "$1" in
    1.as-cast)    echo "restart.glass" ;;
    2.anneal_500) echo "restart.anneal_500K" ;;
    3.anneal_700) echo "restart.anneal_700K" ;;
    *)            echo "" ;;
  esac
}

copy_and_submit_one() {
  local leaf_dir="$1"    # directory that contains lammps_script
  [[ -f "$leaf_dir/lammps_script" ]] || return

  local state_dir rst
  state_dir="$(basename "$leaf_dir")"
  rst="$(restart_for_state "$state_dir")"
  if [[ -z "$rst" ]]; then
    echo "[WARN] Unknown state dir '$state_dir' at '$leaf_dir' -> skipping"
    return
  fi

  mkdir -p "$leaf_dir"
  cp -f "$SRC_DIR/$rst" "$leaf_dir/"   || { echo "[ERR ] Copy failed: $rst -> $leaf_dir"; return; }
  cp -f "$SRC_DIR/slurm_run.sh" "$leaf_dir/" || { echo "[ERR ] Copy failed: slurm_run.sh -> $leaf_dir"; return; }

  ( cd "$leaf_dir" && \
    echo "[SUBM] $(pwd) : using $rst" && \
    jid_out=$(sbatch slurm_run.sh 2>&1) && \
    echo "[SUBM] $jid_out" ) || echo "[ERR ] Submission failed in $leaf_dir"
}

process_tree() {
  local root="$1"
  [[ -d "$root" ]] || return
  while IFS= read -r -d '' script; do
    copy_and_submit_one "$(dirname "$script")"
  done < <(find "$root" -type f -name "lammps_script" -print0 2>/dev/null)
}

main() {
  echo "[INFO] Copy & submit under: ${TARGET_ROOTS[*]}"
  for root in "${TARGET_ROOTS[@]}"; do
    process_tree "$root"
  done
  echo "[DONE] All submissions attempted."
}

main "$@"

