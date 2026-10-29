# Operations: archive, device transfer and disk budget

QStarS moved to an always-on Apple Silicon Mac mini on 2026-10-29 so the closure campaign can run without
pauses. The device has 256 GB; the project (repository + `~/qsrd-runtime`) is kept at **100 GB or less,
including growth**. This page records how results are archived, what is copied to the device, and how
long runs stay inside the budget.

## Where results live

| Place | What | Size (2026-10-28) |
|---|---|---|
| GitHub (`main`) | code, tests, scripts, CMFGEN patches, reports, `Reports/evidence`, benchmark plans | ~0.9 GB |
| `results/` on every machine | the working set: `nk167_trap/fixedT_b35` (stage-1 run, without its set-aside branches), `nk167_trap/exe*`, launch scripts, `b36_tests`, `tprobe_it01`, `nk167_adv/transfer`, `cmfgen_ali_lambda167`, `cmfgen_ali_lambda84_long` | ~11 GB |
| `results/_archive/record/` + `manifests/` + `INDEX.json` (every machine, GitHub release `results-archive-2026-10-29`) | the record of every packed run: inputs, logs, run manifests, metrics, output populations; one manifest line per file with its sha256 | 0.15 + 0.03 GB |
| `results/_archive/full/` (the MacBook; copy to an external disk) | every packed run byte-exact: 321 units, 305.8 GB raw | 113.8 GB |

The device receives 41.4 GB: repository and working set 19.5 GB, MESA 21.1 GB, MESA SDK 0.7 GB,
Claude's project notes 0.15 GB (`transfer_to_device.sh local --dry-run`, 2026-10-29).

Most of the raw volume of a CMFGEN pass is native scratch (`BCSCRATCH*`, `DSCRATCH*`, ~30 % of all
bytes), QStarS operator and spectral traces (compress 7-9x) and radiation-field dumps
(`JH_AT_CURRENT_TIME`, `EDDFACTOR`, `QSTARS_RJ.bin`: double precision, they compress by only ~7 %).
A lossless archive is therefore about a third of the raw size, not a few GB; the record is small because
it leaves those families out.

## Packing runs: `scripts/archive_results.py`

```sh
.venv/bin/python scripts/archive_results.py plan            # units still raw, and the working set
.venv/bin/python scripts/archive_results.py pack --all --delete
.venv/bin/python scripts/archive_results.py status
.venv/bin/python scripts/archive_results.py verify --all     # decompress and re-hash every archive
.venv/bin/python scripts/archive_results.py restore nk167_hot          # full run back into results/
.venv/bin/python scripts/archive_results.py restore nk167_hot --record # the record only
```

A unit is a run directory below `results/` (the working set is never a unit; directories set aside inside
it, such as `fixedT_b35/branch_*`, are). Each unit gives `full/<name>.tar.zst`, `record/<name>.tar.zst`
and `manifests/<name>.jsonl.zst` (`<name>` = path with `/` → `__`). The source is removed only after both
archives were decompressed and every member matched the manifest: type, size, permissions, mtime, link
target and sha256. Restoring by hand needs `zstd` and the window flag:

```sh
zstd -dc --long=30 results/_archive/full/<name>.tar.zst | tar -xpf - -C results
```

The record rule (`RECORD_RULE = 1`) excludes native scratch, QStarS traces, radiation-field and operator
dumps, restart copies of the source populations (`*_IN`), `SCRTEMP`, Fortran unit dumps, rate
diagnostics, `.bin/.npy/.npz` arrays and executables.

## Transfer to the device

Prerequisites on the Mac mini, done by its owner (they need the account password or change system settings):

1. Sign in as `abraiymbek`; create `~/Projects`.
2. System Settings → General → Sharing → **Remote Login** on.
3. System Settings → Energy → prevent automatic sleeping; start up after a power failure.
4. `xcode-select --install`; Homebrew, then `brew install zstd uv gh`.
5. On the MacBook: `ssh-keygen -t ed25519` (if `~/.ssh/id_ed25519` does not exist), then
   `ssh-copy-id abraiymbek@<mac-mini>.local` (asks for the Mac mini password once).
6. On the Mac mini: `gh auth login` (pushing), Claude Code signed in.

Then, from the MacBook:

```sh
scripts/transfer_to_device.sh local --dry-run                     # sizes, no device needed
scripts/transfer_to_device.sh abraiymbek@<mac-mini>.local --dry-run
scripts/transfer_to_device.sh abraiymbek@<mac-mini>.local --verify
```

It refuses a remote with another user, a non-arm64 CPU, no `~/Projects` or less than 60 GB free, copies
the working set to the same absolute paths (every symlink, `history.json`, launch script, CMFGEN Makefile
and the executables' `@rpath` hard-code `/Users/abraiymbek/...`), optionally compares every file by
checksum, and runs `scripts/setup_device.sh` on the device: toolchain and runtime checks, `.venv` from
`requirements-lock.txt` (uv, Python 3.12), the test suite and the disk budget. The MESA run directories
of September in `~/qsrd-runtime` (`boundary-*`, `quest-*`) and CMFGEN builds other than 17, 22 and 31-36
stay on the MacBook.

## Long runs inside the budget: `scripts/retain_runs.py`

The stage-1 driver wrote ~1.7 GB per Newton iterate (about one per hour), roughly half of it scratch.
Without retention the device fills in days. Between or during driver runs:

```sh
.venv/bin/python scripts/retain_runs.py results/nk167_trap/fixedT_b35 --protect it04_trial_l1 --dry-run
.venv/bin/python scripts/retain_runs.py results/nk167_trap/fixedT_b35 --protect it04_trial_l1 \
    --offload /Volumes/<archive-disk>/QStarS_archive     # or abraiymbek@<macbook>.local:<path>
```

It packs pass directories the driver no longer reads (everything except the last `--keep-iterates`
iterates, the first base (the `fixed_metric.py` reference), the bases of the kept history entries and the
next base with their `cols_`, `_from_`, `_tp`, `_ali` directories, `it00_base`, `start_base`, `--protect`
names, and directories written in the last `--idle-min` minutes). `--offload` moves each full archive
to an archive disk; without it the full archives stay in `results/_archive/full` and count against the
budget. `--slim` (lossy, recorded in `QSTARS_SLIM.json`) runs `slim_pass.py` first. Above `--budget-gb`
it exits with status 3 and writes `RETENTION_OVER_BUDGET` in the run directory: start no new iterate until
space is made. Actions are logged in `RETENTION.log`.
