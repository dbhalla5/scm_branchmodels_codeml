#!/usr/bin/env python3

"""
Prepare CODEML branch-site analyses for 50 CYP4V2 migration SCMs.

This script ONLY prepares the CODEML run directories and control files.
It does not execute CODEML.

For each SCM, three analyses are prepared:

    1. ALT model, foreground omega starting at 1
    2. ALT model, foreground omega starting at 2
    3. NULL model, foreground omega fixed at 1

The higher-likelihood ALT solution is selected after CODEML execution.

Input:
    Migration/trees/migration_SCM_###_branchsite.nh

Output:
    Migration/runs/SCM_###/alt/
    Migration/runs/SCM_###/omega_2/
    Migration/runs/SCM_###/null/
"""

from pathlib import Path
import re
import sys


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

BASE = Path(
    "/path/to/SCM/branch_site_model/CYP4V2/Migration"    ############# Change this accordingly. This is for "CYP4V2".
)

TREE_DIR = BASE / "trees"
RUN_DIR = BASE / "runs"

PHY = Path(
    "/path/to/SCM/sub_trees/CYP4V2.phy"                  ############# Change this accordingly. This is for "CYP4V2".
)


# ---------------------------------------------------------------------
# CODEML settings
# ---------------------------------------------------------------------

# These settings correspond to the validated CYP4B1 branch-site setup.

COMMON = {
    "noisy": 3,
    "verbose": 1,
    "runmode": 0,
    "seqtype": 1,
    "CodonFreq": 2,
    "model": 2,
    "NSsites": 2,
    "fix_kappa": 0,
    "kappa": 2,
    "fix_blength": 2,
    "cleandata": 1,
}


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def write_ctl(path, seqfile, treefile, outfile, settings):
    """Write one CODEML control file."""

    lines = [
        f"      seqfile = {seqfile}",
        f"      treefile = {treefile}",
        f"      outfile = {outfile}",
    ]

    for key, value in settings.items():
        lines.append(f"      {key} = {value}")

    path.write_text("\n".join(lines) + "\n")


def find_scm_trees():
    """Find and validate all 50 CYP4V2 branch-site trees."""               ##### change this accordingly. This is for  "CYP4V2".

    trees = []

    for i in range(1, 51):
        tree = TREE_DIR / f"migration_SCM_{i:03d}_branchsite.nh"

        if not tree.is_file():
            raise FileNotFoundError(
                f"Missing SCM tree:\n  {tree}"
            )

        trees.append(tree)

    return trees


# ---------------------------------------------------------------------
# Main preparation
# ---------------------------------------------------------------------

def main():

    print("[INFO] Preparing CYP4V2 branch-site CODEML runs")
    print(f"[INFO] Base directory: {BASE}")
    print(f"[INFO] Tree directory: {TREE_DIR}")
    print(f"[INFO] Run directory: {RUN_DIR}")
    print()

    # ---- Check input alignment ----

    if not PHY.is_file():
        raise FileNotFoundError(
            f"CYP4V2 alignment not found:\n  {PHY}"
        )

    # ---- Check tree directory ----

    if not TREE_DIR.is_dir():
        raise FileNotFoundError(
            f"Tree directory not found:\n  {TREE_DIR}"
        )

    RUN_DIR.mkdir(parents=True, exist_ok=True)

    # ---- Locate all 50 trees ----

    trees = find_scm_trees()

    print(f"[INFO] Found {len(trees)} SCM trees")
    print()

    # ---- Prepare each SCM ----

    for tree in trees:

        match = re.search(r"SCM_(\d{3})_branchsite\.nh$", tree.name)

        if match is None:
            raise ValueError(
                f"Unexpected tree filename: {tree.name}"
            )

        scm = match.group(1)

        scm_dir = RUN_DIR / f"SCM_{scm}"

        alt_dir = scm_dir / "alt"
        omega2_dir = scm_dir / "omega_2"
        null_dir = scm_dir / "null"

        alt_dir.mkdir(parents=True, exist_ok=True)
        omega2_dir.mkdir(parents=True, exist_ok=True)
        null_dir.mkdir(parents=True, exist_ok=True)

        # =============================================================
        # ALT: omega starting value = 1
        # =============================================================

        alt_settings = COMMON.copy()
        alt_settings.update({
            "fix_omega": 0,
            "omega": 1,
        })

        write_ctl(
            alt_dir / "alt.ctl",
            PHY,
            tree,
            alt_dir / f"migration_SCM_{scm}_branchsite_alt.txt",
            alt_settings,
        )

        # =============================================================
        # ALT: omega starting value = 2
        # =============================================================

        omega2_settings = COMMON.copy()
        omega2_settings.update({
            "fix_omega": 0,
            "omega": 2,
        })

        write_ctl(
            omega2_dir / "omega_2.ctl",
            PHY,
            tree,
            omega2_dir / f"migration_SCM_{scm}_branchsite_omega2.txt",
            omega2_settings,
        )

        # =============================================================
        # NULL: omega fixed at 1
        # =============================================================

        null_settings = COMMON.copy()
        null_settings.update({
            "fix_omega": 1,
            "omega": 1,
        })

        write_ctl(
            null_dir / "null.ctl",
            PHY,
            tree,
            null_dir / f"migration_SCM_{scm}_branchsite_null.txt",
            null_settings,
        )

        print(f"[OK] SCM_{scm}")

    # ---- Final summary ----

    print()
    print("[INFO] Preparation completed successfully.")
    print(f"[INFO] SCMs prepared: {len(trees)}")
    print(f"[INFO] Run directory: {RUN_DIR}")
    print()
    print("[INFO] No CODEML analyses were executed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)
