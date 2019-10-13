#!/usr/bin/env bash
set -euo pipefail

# Manually run this script whenever pyproject.toml changes.
#
# This is a work-around for the fact that pip cannot make editable installs 
# from `pyproject.toml`.  For some reason, some web services (namely RTD) seem 
# to require editable installs, so I need to have a `setup.py` file (for the 
# time being).  This script uses flit to generate a `setup.py` file, basically 
# by making a source distribution and extracting the `setup.py` file from it.

tar=$(flit build 2>&1 >/dev/null | grep -oP 'Built sdist: \Kdist/.*.tar.gz')
dir=$(basename ${tar%.tar.gz})

tar -xzf $tar $dir/setup.py
mv $dir/setup.py setup.py
rm -rf $dir

