#!/usr/bin/env python3
"""
LND pipeline runner.

Runs the SageMath pipeline in pipeline/LND_ready.sage without Jupyter.
Sage itself is installed by setup.sh into ./.sage (see README).

Examples:
    python3 run.py --m 5 --r 3 --max-deg 10      # smoke test
    python3 run.py --m 8 --r 5 --max-deg 60
    python3 run.py --all
    python3 run.py --list
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
READY = ROOT / "pipeline" / "LND_ready.sage"
CASES = ROOT / "cases.json"
RUNS = ROOT / "runs"
DEFAULT_Q = ["X"]


def find_sage(explicit: str | None) -> str:
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    env = os.environ.get("SAGE_BIN")
    if env:
        candidates.append(Path(env))
    candidates.append(ROOT / ".sage" / "bin" / "sage")
    which = shutil.which("sage")
    if which:
        candidates.append(Path(which))

    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return str(cand)
    sys.exit(
        "SageMath not found.\n"
        "  Run ./setup.sh first,\n"
        "  or pass --sage /path/to/sage,\n"
        "  or set SAGE_BIN=/path/to/sage."
    )


def build_tail(m: int, r: int, max_deg: int, q: list[str]) -> str:
    q_repr = repr(q)
    return f"""
STOP_AFTER_FILTER = False
print("### RUN m={m} r={r} max_deg={max_deg} Q_list={q_repr} ###", flush=True)
LIST1, LIST2 = run_pipeline({m}, {r}, {max_deg}, Q_list={q_repr})
print("LIST1_LEN", len(LIST1), flush=True)
for _f in LIST1:
    print("LIST1_POLY", _f, flush=True)
print("LIST2_LEN", len(LIST2), flush=True)
for _f in LIST2:
    print("LIST2_POLY", _f, flush=True)
print("### RUN DONE ###", flush=True)
"""


def run_case(sage: str, m: int, r: int, max_deg: int, q: list[str]) -> dict:
    if not READY.is_file():
        sys.exit(f"Missing pipeline file: {READY}")

    RUNS.mkdir(exist_ok=True)
    tag = f"{m}_{r}_{max_deg}"
    header = READY.read_text(encoding="utf-8")
    script = RUNS / f"_run_{tag}.sage"
    script.write_text(header + build_tail(m, r, max_deg, q), encoding="utf-8")

    log_path = RUNS / f"LND_{tag}.log"
    json_path = RUNS / f"LND_{tag}.json"

    cmd = [sage, str(script)]
    print("\n" + "=" * 70, flush=True)
    print(f"[run] D_({m},{r})  max_deg={max_deg}  Q_list={q}", flush=True)
    print(f"[run] {sage} {script}", flush=True)
    print(f"[run] log: {log_path}", flush=True)
    print("=" * 70, flush=True)

    t0 = time.time()
    lines: list[str] = []
    proc = subprocess.Popen(
        cmd,
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    try:
        with log_path.open("w", encoding="utf-8") as log:
            assert proc.stdout is not None
            for raw in proc.stdout:
                line = raw.decode("utf-8", "replace")
                sys.stdout.write(line)
                sys.stdout.flush()
                log.write(line)
                log.flush()
                lines.append(line.rstrip("\r\n"))
    except KeyboardInterrupt:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        print("\n[run] interrupted by user", flush=True)
        raise

    rc = proc.wait()
    wall = round(time.time() - t0, 1)

    list1 = [l.split("LIST1_POLY ", 1)[1] for l in lines if l.startswith("LIST1_POLY ")]
    list2 = [l.split("LIST2_POLY ", 1)[1] for l in lines if l.startswith("LIST2_POLY ")]

    result = {
        "m": m,
        "r": r,
        "max_deg": max_deg,
        "Q_list": q,
        "exit": rc,
        "wall_sec": wall,
        "LIST1": list1,
        "LIST2": list2,
        "log": str(log_path),
    }
    json_path.write_text(json.dumps(result, indent=1), encoding="utf-8")

    print("\n" + "-" * 70, flush=True)
    print(f"[done] exit={rc}  wall={wall}s  LIST1={len(list1)}  LIST2={len(list2)}", flush=True)
    print(f"[done] saved {json_path}", flush=True)
    print("-" * 70, flush=True)
    return result


def load_cases() -> list[dict]:
    if not CASES.is_file():
        sys.exit(f"Missing cases file: {CASES}")
    data = json.loads(CASES.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        sys.exit("cases.json must be a JSON list")
    return data


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the LND SageMath pipeline.")
    parser.add_argument("--m", type=int, help="m parameter of D_(m,r)")
    parser.add_argument("--r", type=int, help="r parameter of D_(m,r)")
    parser.add_argument("--max-deg", type=int, dest="max_deg", help="maximum degree")
    parser.add_argument("--q", default="X", help="comma-separated Q_list, default: X")
    parser.add_argument("--sage", help="path to the sage executable")
    parser.add_argument("--all", action="store_true", help="run every case in cases.json")
    parser.add_argument("--list", action="store_true", help="list cases and exit")
    parser.add_argument("--stop-on-error", action="store_true",
                        help="stop --all at the first failing case")
    args = parser.parse_args()

    if args.list:
        for i, case in enumerate(load_cases(), 1):
            print(f"{i:2d}. m={case['m']} r={case['r']} max_deg={case['max_deg']} "
                  f"q={case.get('q', DEFAULT_Q)}")
        return 0

    sage = find_sage(args.sage)
    print(f"[sage] {subprocess.check_output([sage, '--version'], text=True).strip()}", flush=True)

    if args.all:
        cases = load_cases()
    else:
        if args.m is None or args.r is None or args.max_deg is None:
            parser.error("--m, --r and --max-deg are required unless --all or --list is used")
        q = [x.strip() for x in args.q.split(",") if x.strip()]
        cases = [{"m": args.m, "r": args.r, "max_deg": args.max_deg, "q": q}]

    failures = 0
    for case in cases:
        q = case.get("q", DEFAULT_Q)
        result = run_case(sage, case["m"], case["r"], case["max_deg"], q)
        if result["exit"] != 0:
            failures += 1
            if args.stop_on_error:
                break

    if failures:
        print(f"\n[failed] {failures} case(s) exited non-zero", flush=True)
        return 1
    print("\n[ok] all requested cases finished", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
