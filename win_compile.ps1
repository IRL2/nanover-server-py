# Param statement must be first non-comment, non-blank line in the script
param(
    [switch][alias("e")]$edit = $false,
    [switch][alias("u")]$user = $false
)
    
function announce {
    Write-Host $args[0] -ForegroundColor Green

}

$edit_option = ""
$user_option = "" 

if ($edit)
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

announce "Compiling proto files to python"
python ./python-libraries/narupa-core/setup.py compile_proto

announce "Installing the python packages"
python -m pip install ${edit_option} ${user_option} ./python-libraries/narupa-core/

Get-ChildItem -Directory python-libraries/narupa-* | ForEach-Object { 
    Write-Host "$($_.FullName)"
    pip install ${edit_option} ${user_option} ""$($_.FullName)""
 }

try
{
    python -c "import simtk" 
}
catch 
{
    announce "OpenMM appears to not be installed."
    announce "See <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>."
}

announce "Compiling proto files to C#"
dotnet build --configuration Release csharp-libraries/Narupa.Protocol