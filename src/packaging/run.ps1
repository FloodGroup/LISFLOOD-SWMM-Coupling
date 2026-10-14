# Launch helper (Windows)
# Usage: place lisflood.exe beside this script, then:
#   .\run.ps1 path\to\case.par
param([Parameter(Mandatory=$true)][string]$ParFile)
$exe = Join-Path $PSScriptRoot 'lisflood.exe'
if (-not (Test-Path $exe)) { throw "lisflood.exe not found in $PSScriptRoot" }
& $exe $ParFile
