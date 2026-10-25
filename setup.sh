#!/usr/bin/env bash
# Install a self-contained SageMath into ./.sage (no root, no conda required).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="$ROOT/.sage"
MAMBA_DIR="$ROOT/.micromamba"
SAGE_VERSION="${SAGE_VERSION:-10.9}"

say() { printf '\n[setup] %s\n' "$*"; }
die() { printf '\n[setup] ERROR: %s\n' "$*" >&2; exit 1; }

case "$(uname -s)" in
  Linux)  OS=linux ;;
  Darwin) OS=osx ;;
  *) die "Unsupported OS: $(uname -s). On Windows use WSL2 (see README)." ;;
esac

case "$(uname -m)" in
  x86_64|amd64) ARCH=64 ;;
  aarch64|arm64) ARCH=arm64 ;;
  *) die "Unsupported architecture: $(uname -m)" ;;
esac

if [ "$OS" = "linux" ]; then
  if [ "$ARCH" = "64" ]; then PLATFORM=linux-64; else PLATFORM=linux-aarch64; fi
else
  if [ "$ARCH" = "64" ]; then PLATFORM=osx-64; else PLATFORM=osx-arm64; fi
fi

command -v tar >/dev/null || die "tar is required"
if ! command -v curl >/dev/null && ! command -v wget >/dev/null; then
  die "curl or wget is required"
fi

# Rough disk check: need ~6 GB free.
FREE_KB="$(df -Pk "$ROOT" | awk 'NR==2 {print $4}')"
if [ -n "${FREE_KB:-}" ] && [ "$FREE_KB" -lt 6000000 ]; then
  say "WARNING: only $((FREE_KB / 1024)) MB free in $ROOT; Sage needs ~5-6 GB."
fi

if [ ! -x "$MAMBA_DIR/bin/micromamba" ]; then
  say "Downloading micromamba ($PLATFORM)"
  mkdir -p "$MAMBA_DIR/bin"
  TMP="$(mktemp -d)"
  URL="https://micro.mamba.pm/api/micromamba/$PLATFORM/latest"
  if command -v curl >/dev/null; then
    curl -fL "$URL" -o "$TMP/micromamba.tar.bz2"
  else
    wget -q "$URL" -O "$TMP/micromamba.tar.bz2"
  fi
  tar -xjf "$TMP/micromamba.tar.bz2" -C "$MAMBA_DIR" bin/micromamba
  rm -rf "$TMP"
else
  say "micromamba already present"
fi

export MAMBA_ROOT_PREFIX="$ROOT/.mamba"

if [ ! -x "$ENV_DIR/bin/sage" ]; then
  say "Creating Sage $SAGE_VERSION in $ENV_DIR (5-15 min, ~4 GB)"
  "$MAMBA_DIR/bin/micromamba" create -y \
    -p "$ENV_DIR" \
    -c conda-forge --override-channels \
    "sage=$SAGE_VERSION"
else
  say "Sage environment already exists"
fi

say "Sage version:"
"$ENV_DIR/bin/sage" --version

say "Setup complete."
say "Smoke test:  python3 run.py --m 5 --r 3 --max-deg 10"
say "All cases:   python3 run.py --all"
