# Fresh verified spatial repeat

This archive contains the complete ten-year 256/512/1024/2048-zone
surface-refined backward-Euler hot-envelope campaign, run with finalized,
unchanged package sources. It does **not** validate general quasi-star v1.

- [Campaign and internal-profile norms](20260928T033503.024409Z_spatial_validation/campaign.json)
- [Bitwise comparison with the earlier campaign](20260928T033503.024409Z_spatial_validation/previous_campaign_comparison.json)
- [Convergence figure](20260928T033503.024409Z_spatial_validation/spatial_validation.png)
- [Original-byte archive manifest](archive_manifest.json)
- [Archive verification](archive_verification.json)

Every run reaches ten years and reports `source_snapshot_verified: true`.
All state arrays, time, BH mass, next timestep, step count and complete ledgers
match the earlier campaign bitwise. Its profile comparisons are therefore
confirmed without relying on snapshots taken during active code edits.

The 52 manifest files pass length and SHA-256 checks after decompressing marked
logs. Source snapshots, histories, ledgers, checkpoints, resolved configurations
and first/final profiles are retained. Intermediate profiles remain in the
original local `results/` directory. Original source paths inside reports remain
unchanged for provenance; use the corresponding run suffix within this archive
when reading archived child files.

Reproduce from the repository root in the documented Python environment:

```sh
python scripts/validate_spatial.py --zones 256 512 1024 2048 --end-time 10 --max-dt .02 --mesh log_radius_surface --mechanics backward_euler --max-wall-seconds 900
```

The independent continuum initialization, final internal-profile norms,
boundary mass coverage and separate timestep/remapping/method comparisons are
explained in the [validation report](../../../docs/validation.md).
