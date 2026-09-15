#!/usr/bin/env python3
"""Fetch the published pycwr 1.0.8 source without using the current checkout."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request

VERSION = "1.0.8"
SDIST_URL = (
    "https://files.pythonhosted.org/packages/5a/b2/"
    "9e23cb52dd32a87cb7ae255d28bfa089186ce7ac8d07651a6c85883dc5b7/"
    "pycwr-1.0.8.tar.gz"
)
SDIST_SHA256 = "b6d3d7b872dbe479f8558f9f885bd4da6e7b7458a557e24a4abd0010c896e93f"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--archive", type=Path, help="Use an already downloaded archive, still verifying its hash.")
    parser.add_argument("--extract", action="store_true", help="Extract the verified source for an out-of-tree cibuildwheel build.")
    args = parser.parse_args()
    args.dest.mkdir(parents=True, exist_ok=True)
    archive = args.dest / f"pycwr-{VERSION}.tar.gz"
    if args.archive:
        data = args.archive.read_bytes()
    else:
        request = urllib.request.Request(SDIST_URL, headers={"User-Agent": "pycwr-classroom-wheel-builder"})
        with urllib.request.urlopen(request, timeout=120) as response:
            data = response.read()
    actual = hashlib.sha256(data).hexdigest()
    if actual != SDIST_SHA256:
        raise SystemExit(f"Source checksum mismatch: expected {SDIST_SHA256}, received {actual}")
    archive.write_bytes(data)
    # The release source is unchanged; the smoke test remains outside this tree.
    with tarfile.open(archive, "r:gz") as tar:
        member = tar.getmember(f"pycwr-{VERSION}/PKG-INFO")
        handle = tar.extractfile(member)
        if handle is None:
            raise SystemExit("The release archive does not contain PKG-INFO")
        metadata = handle.read().decode("utf-8")
        if f"\nVersion: {VERSION}\n" not in metadata or "\nName: pycwr\n" not in metadata:
            raise SystemExit("Unexpected release metadata")
        if args.extract:
            source_root = args.dest / f"pycwr-{VERSION}"
            if source_root.exists():
                raise SystemExit(f"Refusing to mix release source with an existing directory: {source_root}")
            for item in tar.getmembers():
                parts = Path(item.name).parts
                if not parts or parts[0] != f"pycwr-{VERSION}" or ".." in parts or Path(item.name).is_absolute():
                    raise SystemExit(f"Unexpected archive path: {item.name}")
                target = args.dest.joinpath(*parts)
                if item.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                elif item.isfile():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    stream = tar.extractfile(item)
                    if stream is None:
                        raise SystemExit(f"Cannot read archive member: {item.name}")
                    target.write_bytes(stream.read())
                else:
                    raise SystemExit(f"Unsupported archive member type: {item.name}")
    record = {"name": "pycwr", "version": VERSION, "url": SDIST_URL, "sha256": actual,
              "bytes": len(data), "source_modified": False}
    (args.dest / "release-provenance.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(archive.resolve())


if __name__ == "__main__":
    main()
