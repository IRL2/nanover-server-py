# Param statement must be first non-comment, non-blank line in the script
param(
    [switch][alias("n")]$noedit = $false,
    [switch][alias("u")]$user = $false
)
    
function announce {
    Write-Host $args[0] -ForegroundColor Green

}

$edit_option = ""
$user_option = "" 

if ($noedit)
{
    announce "Installing narupa in non-edit mode."
}
else
{
    $edit_option = "-e"
    Announce "Installing narupa-protocol in edit mode."
}

if ($user) 
{
    $user_option = "--user"
    Announce "Installing requirements with pip for the user only."
}

announce "Installing python requirements"
python -m pip install -r ./python-libraries/narupa-core/requirements.txt ${user_option}

announce "Installing prototypes requirements"
python -m pip install -r ./python-libraries/prototypes/requirements.txt ${user_option}

announce "Installing python test requirements"
python -m pip install -r ./python-libraries/requirements.test ${user_option}

announce "Compiling proto files to python"
python ./python-libraries/narupa-core/setup.py compile_proto

announce "Installing the python packages"
python -m pip install ${edit_option} ${user_option} ./python-libraries/narupa-core/

Get-ChildItem -Directory python-libraries/narupa-* | ForEach-Object {
    if (Test-Path -Path "$($_.FullName)/setup.py") {
        Write-Host "$($_.FullName)"
        pip install ${edit_option} ${user_option} ""$($_.FullName)""
    }
 }

python -c "import simtk"
if ($LASTEXITCODE -ne 0)
{
    announce "OpenMM appears to not be installed."
    announce "See <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>."
}

python -c "import mpi4py"
if ($LASTEXITCODE -ne 0)
{
    announce "Cannot load mpi4py. Do you have Microsoft MPI installed?"
    announce "See https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi"
}

announce "Compiling proto files to C#"
dotnet build --configuration Release csharp-libraries/Narupa.Protocol
dotnet publish --configuration Release csharp-libraries/Narupa.Protocol