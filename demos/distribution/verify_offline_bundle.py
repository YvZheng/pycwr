#!/usr/bin/env python3
"""Install a classroom wheelhouse offline into a fresh environment and test it."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import venv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    wheelhouse = bundle / "wheels"
    requirements = bundle / "requirements-classroom.txt"
    wheels = sorted(wheelhouse.glob("*.whl"))
    assert wheels and any(p.name.startswith("pycwr-1.0.8-") for p in wheels)
    assert sys.version_info[:2] == (3, 11)
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pycwr-offline-verify-") as temp:
        prefix = Path(temp) / "venv"
        venv.EnvBuilder(with_pip=True, system_site_packages=False).create(prefix)
        python = prefix / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
        run_env = {key: value for key, value in os.environ.items() if not key.startswith("PIP_")}
        # Ignore user package-index settings; all installation candidates must be wheels here.
        run_env["PIP_CONFIG_FILE"] = os.devnull
        run_env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
        run_env.pop("PYTHONPATH", None)
        run_env.pop("PYTHONHOME", None)
        commands = [
            [str(python), "-m", "pip", "install", "-r", "requirements-offline.txt"],
            [str(python), "-m", "pip", "check"],
            [str(python), str(Path(__file__).with_name("smoke_wheel.py").resolve()),
             "--output", str((args.output / "wheel-smoke.json").resolve())],
            [str(python), "-m", "jupyterlab", "--version"],
        ]
        for index, command in enumerate(commands):
            subprocess.run(command, cwd=bundle if index == 0 else temp, env=run_env, check=True)
    record = {"status": "passed", "python": sys.version, "platform": platform.platform(),
              "machine": platform.machine(), "network_used_for_install": False,
              "system_site_packages": False, "source_build_allowed": False,
              "requirements_sha256": hashlib.sha256(requirements.read_bytes()).hexdigest(),
              "files": [{"name": p.name, "bytes": p.stat().st_size,
                         "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in wheels]}
    (args.output / "offline-install.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"Verified {len(wheels)} wheels without an index or source compilation")


if __name__ == "__main__":
    main()
