"""hicstuff declares the packages it imports."""

import ast
import importlib.metadata as md
import pathlib
import re
import sys

import pytest

import hicstuff

SRC = pathlib.Path(hicstuff.__file__).parent


def _imported_modules():
    """Top-level names of the modules hicstuff's source imports."""
    names = set()
    for path in SRC.glob("*.py"):
        # Archived experimental code, no longer supported
        if path.name == "_misc.py":
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                names.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                names.add(node.module.split(".")[0])
    return names


def _normalise(name):
    return re.sub(r"[-_.]+", "-", name).lower()


@pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="packages_distributions finds packages without top_level.txt from Python 3.11",
)
def test_imported_packages_are_declared():
    """Every module hicstuff imports from an installed package comes from a
    package it declares, not from one that another dependency installs."""
    declared = {
        _normalise(re.match(r"[A-Za-z0-9._-]+", req).group(0))
        for req in md.requires("hicstuff")
        if "extra ==" not in req
    }
    packages = md.packages_distributions()
    imported = _imported_modules() & set(packages)
    assert {"numpy", "pysam", "cooler"} <= imported
    undeclared = {
        module: packages[module]
        for module in imported
        if module != "hicstuff" and not declared & {_normalise(p) for p in packages[module]}
    }
    assert not undeclared
