# LND pipeline (portable SageMath runner)

Runs the LND / Maubach / van den Essen pipeline on **any machine with nothing
pre-installed** — no Jupyter, no conda, no admin rights (on Linux/WSL).

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

> On WSL, clone **inside the Linux home** (`cd ~`) rather than `/mnt/c/...`.
> It is much faster and avoids permission quirks.

## Quick start (Windows + WSL2)

1. Install [Git for Windows](https://git-scm.com/download/win) and WSL2
   (`wsl --install -d Ubuntu` as Administrator, then reboot).
2. Clone the repo (anywhere is fine, but WSL home is faster).
3. In PowerShell, from the repo folder:

```powershell
powershell -ExecutionPolicy Bypass -File setup.ps1
```

4. Run a case:

```powershell
wsl -e bash -lc "cd '/mnt/c/Users/<you>/lnd-pipeline' && python3 run.py --m 5 --r 3 --max-deg 10"
```

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
- **Micromamba download fails** → check the network/proxy; the URL is
  `https://micro.mamba.pm/api/micromamba/linux-64/latest`.
- **Disk full** → the env needs ~4 GB; remove `.sage` and re-run setup.
- **`OSError: [Errno 28]` or killed process** → out of memory; close other apps
  or run a smaller `max_deg`.
- **Sage 10.9 vs 9.5** → same pipeline; Gröbner output can differ in ordering
  between Singular versions, but the ideals/relations are the same.
- **No Jupyter is needed.** The notebooks in `notebooks/` are included only for
  reference / the Colab route.

## Colab fallback

If the college PC cannot run WSL2 or install ~4 GB, open
`notebooks/LND_colab_full.ipynb` in Google Colab (CPU runtime). It installs the
same Sage 10.9 environment and runs the same cases.
