#!/usr/bin/env python3
"""Prints the RUN steps of a Dockerfile as a bash script, for builds without Docker.

Usage: dockerfile-steps.py tools/Dockerfile > steps.sh
ARG defaults become exported variables, and ENV values are set for later steps.
"""
import re
import sys

lines = open(sys.argv[1]).read().split("\n")
steps, current, args = [], None, []
for line in lines:
    if current is None:
        match = re.match(r"ARG (\w+)=(\S+)", line)
        if match:
            args.append(f"export {match.group(1)}={match.group(2)}")
        if line.startswith("RUN "):
            current = [line[4:]]
            if not line.rstrip().endswith("\\"):
                steps.append(current)
                current = None
    elif not line.lstrip().startswith("#"):
        current.append(line)
        if not line.rstrip().endswith("\\"):
            steps.append(current)
            current = None

print("#!/bin/bash\nset -euo pipefail")
print("export DEBIAN_FRONTEND=noninteractive TZ=UTC LC_ALL=C.UTF-8 LANG=C.UTF-8 BSPM_DISABLE=1 MAMBA_ROOT_PREFIX=/opt/conda")
print("export PATH=/usr/local/bin:/opt/conda/envs/tools/bin:$PATH")
print("\n".join(args))
for step in steps:
    print("\n".join(step))
