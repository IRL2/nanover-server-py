# Param statement must be first non-comment, non-blank line in the script
param(
    [switch][alias("n")]$noedit = $false
)
    
function announce {
    Write-Host $args[0] -ForegroundColor Green
}

$edit_option = ""

if ($noedit)
{
    announce "Installing nanover in non-edit mode."
}
else
{
    $edit_option = "-e"
    Announce "Installing nanover-server-py in edit mode."
}

announce "Installing the python packages"
python -m pip install ${edit_option} "./python-libraries/nanover-server/[dev]" --config-settings editable_mode=compat

python -c "import openmm"
if ($LASTEXITCODE -ne 0)
{
    announce "OpenMM appears to not be installed."
    announce "See <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>."
}