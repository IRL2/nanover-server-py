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

# Ideally, users should not use "pip install --user". This tends to polute
# the python environment and lead to unpleasant surprises where users are
# using old versions of libraries while they think they installed a newer
# version. Though, using "--user" is sometimes the easiest way of installing
# libraries, so it should be an option.
# Running "./compile.sh --user" will pip install with the --user option.

# A developer most likely want to install narupa's python packages in edit
# mode. If not, they can supply the --no-edit argument.
user_option=""
edit_option="-e"
with_dotnet=true
with_python=true
for option in "$@"; do
    if [[ "$option" == "--user" ]]; then
        user_option="--user"
    elif [[ "$option" == "--no-edit" ]]; then
        edit_option=""
    elif [[ "$option" == "--no-dotnet" ]]; then
        with_dotnet=false
    elif [[ "$option" == "--no-python" ]]; then
        with_python=false
    fi
done
# We do not want to use pip with --user if we use -e.
narupa_user_option=${user_option}
if [[ ! -z "${edit_option}" ]]; then
    narupa_user_option=""
fi

if [[ $with_python == true ]]; then
    # mpi4py (required for narupa-lammps) needs MPI to be installed on the system.
    python -c "import mpi4py" 2>&1 > /dev/null || {
        announce "The mpi4py library is required but cannot be found."
        announce "Because it requires a system library, it is not installed by this"
        announce "script; you have to install it yourself."
        announce "If you are using conda, install mpi4py by running"
        announce "conda install -c conda-forge mpi4py"
        exit 1
    }


    announce "Installing python requirements"
    python -m pip install -r ./python-libraries/narupa-core/requirements.txt ${user_option}

    announce "Installing python prototypes requirements"
    python -m pip install -r ./python-libraries/prototypes/requirements.txt ${user_option}

    announce "Installing python test requirements"
    python -m pip install -r ./python-libraries/requirements.test ${user_option}

    announce "Compiling proto files to python"
    python ./python-libraries/narupa-core/setup.py compile_proto

    announce "Installing the python packages"
    python -m pip install ${edit_option} ${narupa_user_option} ./python-libraries/narupa-core/

    for package in python-libraries/narupa-*/; do
        if [[ -f "${package}/setup.py" ]]; then
            python -m pip install ${edit_option} ${narupa_user_option} ${package}
        fi
    done

    python -c "import simtk" 2>&1 > /dev/null || {
        announce "OpenMM is not installed."
        announce "See <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>."
    }
fi

if [[ $with_dotnet == true ]]; then
    announce "Compiling proto files to C#"
    dotnet build --configuration Release csharp-libraries/Narupa.Protocol
    dotnet publish --configuration Release csharp-libraries/Narupa.Protocol
fi
