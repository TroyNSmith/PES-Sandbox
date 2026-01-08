from dataclasses import dataclass
from pathlib import Path

from math import floor
import copy

@dataclass
class ORCA_Parameters:
    # == I/O
    xyz_in: str | Path
    name_out: str = "calc"

    # == Input parameters
    multiplicity: int = 1
    functional: str = ""
    basis: str = ""
    method_keywords: str = ""
    block_inputs: str = ""

    # == Submit parameters
    bash_commands: str = ""

    # == Resource requirements
    processors: int = 8
    max_memory: int = 1000
    lscratch_size: int = 20
    time_limit: str = "04:00:00"


def orca_inputs(
    amchi: str,
    pars: ORCA_Parameters,
    data_dir: str | Path,
):
    """Writes .inp and .sh files for ORCA calculation."""
    pars = copy.deepcopy(pars)

    inp_text = f"""%PAL NPROCS {pars.processors} END
%MaxCore {floor(pars.max_memory * 0.8)}
%base "{pars.name_out}"
! {pars.functional} {pars.basis} {pars.method_keywords}
{pars.block_inputs}
* xyzfile 0 {pars.multiplicity} {pars.xyz_in}
"""

    sh_header = """#!/bin/bash
set -euo pipefail
SUBMIT_DIR=$(pwd)
echo "SUBMITTED: $(date)"

SCRATCH_DIR="/lscratch/${USER}/${RANDOM}"
echo "SCRATCH_DIR: $SCRATCH_DIR"
echo "NODE: $(hostname)"

mkdir -p "$SCRATCH_DIR"
rsync -a "$SUBMIT_DIR/" "$SCRATCH_DIR/"
cd "$SCRATCH_DIR"

cleanup() {
    rsync -a "$SCRATCH_DIR/" "$SUBMIT_DIR/"
    rm -rf "$SCRATCH_DIR"
    cd "$SUBMIT_DIR"
}
trap cleanup EXIT\n
"""

    sh_body = f"""module load ORCA/6.1
$(which orca) {pars.name_out}.inp > "$SUBMIT_DIR/{pars.name_out}.log"
{pars.bash_commands}
"""
    if isinstance(data_dir, str):
        data_dir = Path(data_dir)

    directory = data_dir / amchi
    (directory / f"{pars.name_out}.inp").write_text(inp_text)
    (directory / f"{pars.name_out}.sh").write_text(sh_header + sh_body)

    return directory / f"{pars.name_out}.sh"
