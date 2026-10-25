# LND pipeline (portable SageMath runner)

Runs the LND / Maubach / van den Essen pipeline on **any machine with nothing
pre-installed** — no Jupyter, no conda. Linux/macOS need no admin rights;
Windows needs Administrator **once** to enable WSL2.

`setup.sh` downloads [micromamba](https://micro.mamba.pm) and installs a
self-contained **SageMath 10.9** into `./.sage` (about 4 GB, 5–15 min).
`run.py` then runs your pipeline in a Sage subprocess, streams every line of
output to the terminal, and saves the full log plus `LIST1`/`LIST2` as JSON.

---

## What you need

| | Requirement |
|---|---|
| OS | Linux (x86_64 / aarch64), macOS, or **Windows via WSL2** |
| RAM | 8 GB minimum, 16–20 GB recommended (your 20 GB college PC is good) |
| Disk | ~6 GB free for the Sage environment |
| Software | `git` + `python3` (already present on most Linux/WSL installs) |
| GPU | **not used** — Sage/Singular Gröbner bases are CPU-only |

> **Windows:** SageMath has no native Windows build. You need WSL2 + Ubuntu.
> In an **Administrator** PowerShell: `wsl --install -d Ubuntu`, reboot, then
> run `setup.ps1`. If the college PC blocks WSL, use the Colab notebook instead
> (`notebooks/LND_colab_full.ipynb`).

---

## Quick start (Linux / WSL / macOS)

```bash
git clone https://github.com/vihan17/lnd-pipeline.git
cd lnd-pipeline

bash setup.sh        # installs Sage 10.9 into ./.sage (once)

python3 run.py --m 5 --r 3 --max-deg 10     # smoke test, ~seconds
python3 run.py --m 8 --r 5 --max-deg 60     # your target case
python3 run.py --all                        # all cases in cases.json
```

> **Public repo** — `git clone` needs no GitHub login.
>
> On WSL, clone **inside the Linux home** (`cd ~`) rather than `/mnt/c/...`.
> It is much faster and avoids permission quirks.

## Quick start (Windows + WSL2)

SageMath has no native Windows build, so everything runs inside Ubuntu (WSL2).
You need Administrator rights **once** to enable WSL2.

### 1. Enable WSL2 (Administrator PowerShell, once)

```powershell
wsl --install -d Ubuntu
```

Reboot, then open the **Ubuntu** app and create your Linux username/password.
If `wsl --install` fails, enable virtualization (VT-x/AMD-V) in the BIOS, or
ask the lab admin to enable the "Windows Subsystem for Linux" feature.

### 2. Install everything inside Ubuntu (not PowerShell)

```bash
cd ~
git clone https://github.com/vihan17/lnd-pipeline.git
cd lnd-pipeline

sudo apt update && sudo apt install -y git   # only if git is missing
bash setup.sh                                 # 10-15 min, ~4 GB
python3 run.py --m 5 --r 3 --max-deg 10       # smoke test
```

### 3. Optional: give WSL more RAM

WSL2 defaults to ~50% of RAM (~10 GB on a 20 GB PC). If the Gröbner stage needs
more, create `C:\Users\<you>\.wslconfig`:

```ini
[wsl2]
memory=14GB
processors=8
```

Then run `wsl --shutdown` in PowerShell and reopen Ubuntu.

> Already cloned on the Windows side (`C:\...`)? You can run
> `powershell -ExecutionPolicy Bypass -File setup.ps1`; it calls `setup.sh`
> through WSL. It works, but the Sage env then lives on `/mnt/c`, which is
> noticeably slower than cloning inside `~`.

---

## Commands

```bash
python3 run.py --list                  # show cases from cases.json
python3 run.py --all                   # run every case, in order
python3 run.py --m 8 --r 5 --max-deg 60
python3 run.py --m 5 --r 3 --max-deg 38 --q X
python3 run.py --sage /usr/bin/sage --m 5 --r 3 --max-deg 10
```

Options:

- `--q X,Y` — variables for `Q_list` (default `X`; use `--q ""` for none)
- `--sage PATH` — use an existing Sage instead of `./.sage/bin/sage`
- `--stop-on-error` — stop `--all` at the first non-zero exit

Want the notebook UI? After setup:

```bash
./.sage/bin/sage -n jupyter --no-browser
```

Open the printed URL and pick `notebooks/LND.ipynb` (Sage 10.9 kernel).

## Outputs

Everything lands in `runs/`:

```
runs/LND_8_5_60.log     # complete console output (same as Jupyter)
runs/LND_8_5_60.json    # {m, r, max_deg, Q_list, LIST1, LIST2, exit, wall_sec}
runs/_run_8_5_60.sage   # exact script that was executed
```

`LIST1` / `LIST2` are stored as lists of polynomial strings, exactly as the
notebook prints them.

## Cases in `cases.json`

| m | r | max_deg | note |
|---|---|---------|------|
| 5 | 3 | 38 | lighter real case |
| 8 | 3 | 40 | |
| 8 | 5 | 47 | known ~24 s locally (baseline) |
| 8 | 5 | 55 | scaling probe |
| 8 | 5 | 60 | heaviest |

Run them in this order. If `47` already takes very long on the college PC,
`60` will not finish quickly — the bottleneck is the Gröbner/membership stage,
not the monomial enumeration.

## How it works

1. `pipeline/LND_ready.sage` contains the full pipeline definitions (the same
   code as the notebook cells, with `STOP_AFTER_FILTER` overridden at runtime).
2. `run.py` appends a small tail that calls
   `run_pipeline(m, r, max_deg, Q_list=[...])`, prints `LIST1_POLY` /
   `LIST2_POLY` markers, and runs the result with the Sage binary.
3. Output is tee'd to the terminal and to `runs/*.log`; the markers are parsed
   into JSON.

## Notes / troubleshooting

- **`SageMath not found`** → run `bash setup.sh` first.
- **WSL: `Could not resolve host`** → run `wsl --shutdown` in PowerShell, reopen
  WSL and retry. WSL networking occasionally loses DNS after sleep.
- **Windows without Administrator rights** → WSL2 cannot be enabled; use the
  Colab notebook instead (`notebooks/LND_colab_full.ipynb`).
- **Micromamba download fails** → check the network/proxy; the URL is
  `https://micro.mamba.pm/api/micromamba/linux-64/latest` (7 MB), after which
  conda-forge packages (~1.5 GB download, ~4 GB installed) are fetched.
- **Disk full** → the env needs ~4 GB; remove `.sage` and re-run setup.
- **`OSError: [Errno 28]` or killed process** → out of memory; close other apps
  or run a smaller `max_deg`.
- **Sage 10.9 vs 9.5** → same pipeline; Gröbner output can differ in ordering
  between Singular versions, but the ideals/relations are the same.
- **No Jupyter is needed.** The notebooks in `notebooks/` are included only for
  reference / the Colab route.

## Colab fallback

If the PC cannot run WSL2 (no admin rights) or cannot spare ~4 GB, open the
notebook directly from this public repo:

**https://colab.research.google.com/github/vihan17/lnd-pipeline/blob/main/notebooks/LND_colab_full.ipynb**

Choose a **CPU** runtime, run the **Setup** cell once (5–15 min), then the case
cells. It installs the same Sage 10.9 and runs the same pipeline.
