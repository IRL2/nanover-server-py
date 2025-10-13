param(
    [switch][alias("n")]$noedit = $false
)
    
function announce {
    Write-Host $args[0] -ForegroundColor Green
}

if ($noedit)
{
    announce "Installing nanover in non-edit mode."
    python -m pip install "./python-libraries/nanover-server/[dev]"
}
else
{
    announce "Installing nanover-server-py in edit mode."
    python -m pip install -e "./python-libraries/nanover-server/[dev]" --config-settings editable_mode=compat
}

python -c "import openmm"
if ($LASTEXITCODE -ne 0)
{
    announce "OpenMM appears to not be installed."
    announce "See <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>."
}