#!/usr/bin/env python3
"""Run the wheel test suite without shell glob expansion or source imports."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import unittest

def main():
    from ci_import_smoke import _sanitize_sys_path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    _sanitize_sys_path(root)

    import pycwr

    module_path = Path(pycwr.__file__).resolve()
    if module_path.is_relative_to(root):
        raise RuntimeError(f"Tests must exercise the installed wheel: {module_path}")
    suite = unittest.defaultTestLoader.discover(str(root / "test"), pattern="test_*.py")
    count = suite.countTestCases()
    if count < 293:
        raise RuntimeError(f"Incomplete discovery: expected at least 293 tests, found {count}")
    print(f"Discovered {count} tests against {module_path}", flush=True)
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
