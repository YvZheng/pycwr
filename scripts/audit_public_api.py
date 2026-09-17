#!/usr/bin/env python3
"""Map public Python functions to measured coverage and test contexts.

This reports execution coverage, not a proof that every input is correct.
Run coverage with dynamic_context=test_function and export JSON with contexts.
"""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path


def public_functions(tree, prefix=""):
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            if not node.name.startswith("_"):
                yield from public_functions(node, prefix + node.name + ".")
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if not node.name.startswith("_") or node.name == "__init__":
                yield prefix + node.name, node


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coverage_json", type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    coverage = json.loads(args.coverage_json.read_text())
    entries = []
    for path in sorted((root / "pycwr").rglob("*.py")):
        relative = path.relative_to(root).as_posix()
        stats = coverage["files"].get(relative, coverage["files"].get(str(path), {}))
        executed = set(stats.get("executed_lines", []))
        missing = set(stats.get("missing_lines", []))
        statements = executed | missing
        contexts = stats.get("contexts", {})
        module = relative[:-3].replace("/", ".")
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for symbol, node in public_functions(tree):
            body_lines = set(range(node.body[0].lineno, node.end_lineno + 1))
            # A docstring is executed at definition time and does not prove a call.
            if isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
                body_lines -= set(range(node.body[0].lineno, node.body[0].end_lineno + 1))
            relevant = statements & body_lines
            covered = executed & relevant
            tests = sorted({context for line in covered for context in contexts.get(str(line), []) if context})
            entries.append({
                "symbol": module + "." + symbol,
                "file": relative,
                "line": node.lineno,
                "executed_statements": len(covered),
                "statements": len(relevant),
                "called": bool(covered),
                "missing_lines": sorted(missing & relevant),
                "test_contexts": tests,
            })
    summary = {
        "public_functions": len(entries),
        "called_functions": sum(entry["called"] for entry in entries),
        "uncalled_functions": sum(not entry["called"] for entry in entries),
        "coverage_totals": coverage.get("totals", {}),
        "scope": "Public Python functions/methods and constructors; Cython is verified separately by installed-wheel tests.",
        "limitation": "Execution coverage does not prove correctness for every input or file format.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"summary": summary, "functions": entries}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
