param(
    [string]$SourceDirectory = (Join-Path (Split-Path $PSScriptRoot -Parent) 'src'),
    [string]$OutputDirectory = (Join-Path (Split-Path $PSScriptRoot -Parent) 'build\release'),
    [string]$BuildDirectory = (Join-Path (Split-Path $PSScriptRoot -Parent) 'build\nmake'),
    [string]$VcVars = '<VCVARS64_BAT>',
    [string]$CMake = '<CMAKE_EXE>',
    [string]$CudaRoot = '',
    [int]$CudaArchitecture = 75,
    [string]$MsvcToolset = '14.39',
    [string]$RedistDirectory = ''
)
$ErrorActionPreference = 'Stop'
if ($VcVars -like '<*' -or $CMake -like '<*') {
    throw 'Set -VcVars and -CMake to your local toolchain paths before building.'
}
if ($CudaRoot) {
    $nvcc = Join-Path $CudaRoot 'bin\nvcc.exe'
    if (-not (Test-Path -LiteralPath $nvcc)) { throw "CUDA compiler not found: $nvcc" }
    $env:CUDA_PATH = $CudaRoot
    $env:PATH = (Join-Path $CudaRoot 'bin') + ';' + $env:PATH
}
New-Item -ItemType Directory -Force -Path $BuildDirectory, $OutputDirectory | Out-Null
$envLines = & $env:ComSpec /d /c "call `"$VcVars`" -vcvars_ver=$MsvcToolset >nul && set"
if ($LASTEXITCODE -ne 0) { throw 'MSVC environment initialization failed' }
foreach ($line in $envLines) {
    if ($line -match '^([^=]+)=(.*)$') {
        [Environment]::SetEnvironmentVariable($Matches[1], $Matches[2], 'Process')
    }
}
$configFile = Join-Path $SourceDirectory 'config.default.cmake'
$configureArgs = @(
    '-S', $SourceDirectory,
    '-B', $BuildDirectory,
    '-G', 'NMake Makefiles',
    '-DCMAKE_BUILD_TYPE=Release',
    "-DCMAKE_CUDA_ARCHITECTURES=$CudaArchitecture",
    "-D_CONFIG=$configFile"
)
if ($CudaRoot) {
    $configureArgs += "-DCMAKE_CUDA_COMPILER=$($nvcc.Replace('\', '/'))"
    $configureArgs += "-DCUDAToolkit_ROOT=$($CudaRoot.Replace('\', '/'))"
}
& $CMake @configureArgs
if ($LASTEXITCODE -ne 0) { throw 'CMake configuration failed' }
& $CMake --build $BuildDirectory --target lisflood --config Release
if ($LASTEXITCODE -ne 0) { throw 'Compilation failed' }
foreach ($name in @('lisflood.exe','netcdf.dll','hdf5.dll','hdf5_hl.dll','zlib.dll','msvcr100.dll')) {
    $src = Join-Path $BuildDirectory $name
    if (Test-Path -LiteralPath $src) {
        Copy-Item -LiteralPath $src -Destination $OutputDirectory -Force
    }
}
if ($RedistDirectory -and (Test-Path -LiteralPath $RedistDirectory)) {
    foreach ($name in @('msvcp140.dll','vcruntime140.dll','vcruntime140_1.dll')) {
        $p = Join-Path $RedistDirectory "Microsoft.VC143.CRT\$name"
        if (Test-Path -LiteralPath $p) { Copy-Item $p $OutputDirectory -Force }
    }
    $omp = Join-Path $RedistDirectory 'Microsoft.VC143.OpenMP\vcomp140.dll'
    if (Test-Path -LiteralPath $omp) { Copy-Item $omp $OutputDirectory -Force }
}
Write-Host "Build output: $OutputDirectory"
