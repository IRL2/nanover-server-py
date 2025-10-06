#!/bin/bash

function announce() {
	echo -e "\e[1;34m$1\e[0m"
}

# Stop if something goes wrong to avoid building on errors.
# * "set -e" mean that the script stops if a commands ends with a non-zero
#   exit code
# * "set -u" stops if we refer to an undefined variable
# * "set -o pipefail" considers that a series of piped command fails not only
#   if the last command fails, but if any of them do
set -euo pipefail

# A developer most likely want to install nanover's python packages in edit
# mode. If not, they can supply the --no-edit argument.
edit_option="-e"
with_python=true
for option in "$@"; do
	if [[ "$option" == "--no-edit" ]]; then
		edit_option=""
	fi
done

if [[ $with_python == true ]]; then
	announce "Installing the python packages"
	python -m pip install ${edit_option} ./python-libraries/nanover-server/[dev] --config-settings editable_mode=compat

	python -c "import openmm" 2>&1 >/dev/null || {
		announce "OpenMM is not installed."
		announce "See <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>."
	}
fi
