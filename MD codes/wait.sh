#!/usr/bin/env bash
# Wait for a SLURM job to finish, then copy restarts & submit all leaf jobs.
# Run this from a *composition* directory that has 1.glass/{restart.*,slurm_run.sh}

JOBID="${JOBID:-1574587}"        # override with: JOBID=123456 ./wait_copy_submit.sh
SRC_DIR="1.glass"
REQ_FILES=( "restart.glass" "restart.anneal_500K" "restart.anneal_700K" "slurm_run.sh" )
TARGET_ROOTS=( "2.indent" "3.srjt" "4.creep" )

# --- helper: wait for SLURM job to finish ---
wait_for_job() {
  local jid="$1"
  echo "[INFO] Waiting for SLURM job $jid to finish..."
  # Poll squeue; when it disappears, proceed.
  while true; do
    if [[ -z "$(squeue -h -j "$jid" 2>/dev/null)" ]]; then
      break
    fi
    sleep 10
  done
  # Report final state if sacct is available
  if command -v sacct >/dev/null 2>&1; then
    local st
    st=$(sacct -j "$jid" --format=JobID,State -n -P 2>/dev/null | awk -F'|' 'NR==1{print $2}')
    [[ -n "$st" ]] && echo "[INFO] Job $jid final state: $st"
  else
    echo "[INFO] Job $jid no longer in squeue."
  fi
}

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
  local state_dir base rst

  if [[ ! -f "$leaf_dir/lammps_script" ]]; then
    return
  fi

  state_dir="$(basename "$leaf_dir")"
  rst="$(restart_for_state "$state_dir")"

  if [[ -z "$rst" ]]; then
    echo "[WARN] Unknown state dir '$state_dir' at '$leaf_dir' -> skipping"
    return
  fi

  mkdir -p "$leaf_dir"

  # Copy the appropriate restart + slurm_run.sh
  cp -f "$SRC_DIR/$rst" "$leaf_dir/" || { echo "[ERR ] Copy failed: $rst -> $leaf_dir"; return; }
  cp -f "$SRC_DIR/slurm_run.sh" "$leaf_dir/" || { echo "[ERR ] Copy failed: slurm_run.sh -> $leaf_dir"; return; }

  # Submit
  ( cd "$leaf_dir" && \
    echo "[SUBM] $(pwd) : using $rst" && \
    jid_out=$(sbatch slurm_run.sh 2>&1) && \
    echo "[SUBM] $jid_out" ) || echo "[ERR ] Submission failed in $leaf_dir"
}

process_tree() {
  local root="$1"
  [[ -d "$root" ]] || return
  # Find every directory that has a file named 'lammps_script'
  while IFS= read -r -d '' script; do
    local leaf_dir
    leaf_dir="$(dirname "$script")"
    copy_and_submit_one "$leaf_dir"
  done < <(find "$root" -type f -name "lammps_script" -print0 2>/dev/null)
}

main() {
  wait_for_job "$JOBID"
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
  echo "[INFO] Starting copy & submit under: ${TARGET_ROOTS[*]}"
  for root in "${TARGET_ROOTS[@]}"; do
    process_tree "$root"
  done
  echo "[DONE] All submissions attempted."
}

main "$@"

