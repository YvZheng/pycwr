#!/usr/bin/env python3
"""Smoke-check that an installed pycwr package imports and runs."""

from __future__ import annotations

import argparse
import pkgutil
import sys
from pathlib import Path


EXPECTED_SUBMODULES = [
    "GraphicalInterface",
    "configure",
    "core",
    "draw",
    "interp",
    "io",
    "qc",
    "retrieve",
]


def _sanitize_sys_path(repo_root: Path | None) -> None:
    blocked = set()
    if repo_root is not None:
        blocked.add(str(repo_root.resolve()))
        blocked.add(str((repo_root / "scripts").resolve()))
    blocked.add(str(Path.cwd().resolve()))
    sys.path[:] = [entry for entry in sys.path if entry and str(Path(entry).resolve()) not in blocked]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=None, help="Source-tree root to exclude from sys.path")
    parser.add_argument("--sample", default=None, help="Optional radar sample file for a tiny runtime smoke test")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve() if args.repo_root else None
    _sanitize_sys_path(repo_root)

    import numpy as np
    import pycwr
    import pycwr.GraphicalInterface
    import pycwr.configure
    import pycwr.core
    import pycwr.draw
    import pycwr.interp
    import pycwr.io
    import pycwr.qc
    import pycwr.retrieve
    from pycwr.io import read_auto

    module_file = Path(pycwr.__file__).resolve()
    if repo_root is not None and module_file.is_relative_to(repo_root):
        raise RuntimeError(f"pycwr imported from source tree instead of site-packages: {module_file}")

    discovered = sorted(module.name for module in pkgutil.iter_modules(pycwr.__path__))
    missing = sorted(set(EXPECTED_SUBMODULES) - set(discovered))
    if missing:
        raise RuntimeError(f"installed package is missing top-level submodules: {missing}")

    result = {
        "version": pycwr.__version__,
        "module_file": str(module_file),
        "submodules": discovered,
    }

    if args.sample:
        sample = Path(args.sample).resolve()
        radar = read_auto(str(sample))
        x = np.array([-10_000.0, 0.0, 10_000.0], dtype=np.float64)
        y = np.array([-10_000.0, 0.0, 10_000.0], dtype=np.float64)
        radar.add_product_CR_xy(x, y)
        result["sample_nsweeps"] = int(radar.nsweeps)
        result["has_cr_native"] = "CR_native" in radar.product

    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
