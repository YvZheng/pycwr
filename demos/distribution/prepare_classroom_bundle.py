#!/usr/bin/env python3
"""Put verified-platform wheels and short offline-install requirements together."""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheelhouse", type=Path, required=True)
    parser.add_argument("--requirements", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args()
    wheels = sorted(args.wheelhouse.glob("*.whl"))
    if not wheels or not any(p.name.startswith("pycwr-1.0.8-") for p in wheels):
        raise SystemExit("Missing the pycwr 1.0.8 binary wheel")
    if args.dest.exists():
        raise SystemExit(f"Refusing to mix bundles in an existing directory: {args.dest}")
    target = args.dest / "wheels"
    target.mkdir(parents=True)
    for wheel in wheels:
        shutil.copy2(wheel, target / wheel.name)
    shutil.copy2(args.requirements, args.dest / "requirements-classroom.txt")
    (args.dest / "requirements-offline.txt").write_text(
        "--no-index\n--only-binary=:all:\n--find-links ./wheels\n-r requirements-classroom.txt\n",
        encoding="utf-8",
    )
    print(f"Prepared {len(wheels)} wheels in {args.dest.resolve()}")


if __name__ == "__main__":
    main()
