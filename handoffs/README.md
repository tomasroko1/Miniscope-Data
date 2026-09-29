# CA1 cluster handoff, 29 September 2026

This branch contains a small, self-contained handoff for the R005 registration follow-up and the auditable six-day S–D coactivity export. Start with [`CLUSTER_HANDOFF_2026_09_29/PROMPT_REANUDAR_CLUSTER.md`](CLUSTER_HANDOFF_2026_09_29/PROMPT_REANUDAR_CLUSTER.md).

The `carrillo_poster/` subfolder adds two ready-to-use Stoixeion figures and a separate short explanation of what that method found. Its core-overlap inference comes from the 999-permutation phasewise run; the six-day pairwise-network result elsewhere in the handoff is a different analysis.

To download **only this handoff** into the cluster results directory without switching the existing checkout:

```bash
git -C /mnt/NAS/Tomas/Miniscope-Data fetch origin handoff/r005-ca1-20260929
git -C /mnt/NAS/Tomas/Miniscope-Data archive FETCH_HEAD handoffs/CLUSTER_HANDOFF_2026_09_29 | tar -x -C /mnt/NAS/Tomas/results --strip-components=1
```

The files will be in `/mnt/NAS/Tomas/results/CLUSTER_HANDOFF_2026_09_29/`. Check `SHA256_MANIFEST.csv` before using the copied sources. The `network_10day_reported` directory holds an aggregate report whose detailed four source CSVs still need to be found; `network_6day` has the per-day CSVs and protocol available now.
