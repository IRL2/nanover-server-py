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
    announce "Installing nanover in non-edit mode."
}
else
{
    $edit_option = "-e"
    Announce "Installing nanover-server-py in edit mode."
}

if ($user) 
{
    $user_option = "--user"
    Announce "Installing requirements with pip for the user only."
}

announce "Installing python requirements"
python -m pip install -r ./python-libraries/nanover-server/requirements.txt ${user_option}

announce "Installing prototypes requirements"
python -m pip install -r ./python-libraries/prototypes/requirements.txt ${user_option}

announce "Installing python test requirements"
python -m pip install -r ./python-libraries/requirements.test ${user_option}

announce "Compiling proto files to python"
python ./python-libraries/compile_proto.py --proto-dir=./protocol --python-dir=./python-libraries/nanover-core/src

announce "Installing the python packages"
python -m pip install ${edit_option} ${user_option}  (Convert-Path "./python-libraries/nanover-server/") --config-settings editable_mode=compat

Get-ChildItem -Directory python-libraries/nanover-* | ForEach-Object {
    if (Test-Path -Path "$($_.FullName)/pyproject.toml") {
        Write-Host "$($_.FullName)"
        pip install ${edit_option} ${user_option} ""$($_.FullName)""  --config-settings editable_mode=compat
    }
 }

python -c "import openmm"
if ($LASTEXITCODE -ne 0)
{
    announce "OpenMM appears to not be installed."
    announce "See <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>."
}

announce "Compiling proto files to C#"
dotnet build --configuration Release csharp-libraries/NanoVer.Protocol
dotnet publish --configuration Release csharp-libraries/NanoVer.Protocol