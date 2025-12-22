#!/bin/bash
$PYTHON -m pip install uv
$PYTHON -m uv pip install . --no-deps --ignore-installed --no-cache-dir -vvv