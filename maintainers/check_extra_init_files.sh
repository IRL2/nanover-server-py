#!/bin/bash

# Test for __init__.py files at the root of the nanover python packages. NanoVer,
# on the python side, is a namespace package (see
# <https://packaging.python.org/guides/packaging-namespace-packages/>). For
# this to work, there should NOT be any __init__.py files at the root of the
# subpackages: a __init__.py file is python-libraries/nanover-*/src/nanover will
# make at least this subpackage unable to be imported.

# This script search for these superfluous __init__.py files and exits with a
# return code different from 0 if any is found.

LS_COMMAND="ls python-libraries/nanover-*/src/nanover/__init__.py"

exit_code=0
n_init_files=$(${LS_COMMAND} 2> /dev/null | tee | wc -l)
if [[ ${n_init_files} > 0 ]]; then
    echo "Check the root of the python subpackages for superfluous __init__.py files."
    $LS_COMMAND
    exit_code=1
fi
exit ${exit_code}
