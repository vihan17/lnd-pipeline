# Run the LND pipeline from Windows using WSL2.
# Requires WSL2 + a Linux distro (Ubuntu). This script only orchestrates;
# all installation happens inside WSL.

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    Write-Host "WSL is not installed on this machine." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Open PowerShell as Administrator and run:"
    Write-Host "    wsl --install -d Ubuntu"
    Write-Host ""
    Write-Host "Reboot, create your Linux user, then run this script again."
    exit 1
}

wsl.exe -e true 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "No working WSL2 distribution found." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Open PowerShell as Administrator and run:"
    Write-Host "    wsl --install -d Ubuntu"
    Write-Host ""
    Write-Host "Reboot, create your Linux user, then run this script again."
    exit 1
}

$wslPath = (wsl.exe wslpath -a "$root").Trim()
Write-Host "Repo inside WSL: $wslPath"

if ($wslPath -like "/mnt/*") {
    Write-Host ""
    Write-Host "WARNING: the repo is on a Windows drive (/mnt/...)." -ForegroundColor Yellow
    Write-Host "The Sage environment will be slow there and can hit permission issues."
    Write-Host "Recommended (inside the Ubuntu app):"
    Write-Host "    cd ~"
    Write-Host "    git clone https://github.com/vihan17/lnd-pipeline.git"
    Write-Host "    cd lnd-pipeline && bash setup.sh"
    Write-Host ""
}

Write-Host "Running setup.sh (10-15 min, ~4 GB) ..."

wsl.exe -e bash -lc "cd '$wslPath' && bash setup.sh"
if ($LASTEXITCODE -ne 0) {
    Write-Host "Setup failed (exit $LASTEXITCODE)." -ForegroundColor Red
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host "Setup complete. Smoke test:" -ForegroundColor Green
Write-Host "    wsl -e bash -lc `"cd '$wslPath' && python3 run.py --m 5 --r 3 --max-deg 10`""
Write-Host ""
Write-Host "Full case list:"
Write-Host "    wsl -e bash -lc `"cd '$wslPath' && python3 run.py --all`""
