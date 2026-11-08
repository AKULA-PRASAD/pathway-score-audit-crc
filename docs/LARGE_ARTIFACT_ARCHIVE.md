# Archiving the large Phase 9 artifacts

**Problem.** GitHub rejects any file over 100 MB in a normal Git push, and warns above 50 MB. At B_out = 750 each corrected checkpoint `reports/phase9/corrected_rerun/runs/run<r>.jsonl` is about 113 MB, five of them about 565 MB in total. They are the primary scientific evidence and must be preserved byte for byte.

**Principle.** Nothing is deleted, truncated or altered. The checkpoints stay where the pipeline wrote them. Git holds a SHA-256 manifest that pins every archived byte, and the bytes themselves go to an archive that accepts large files. Every verification step (`verify_sizing_decision`, the post-run verifier, `scripts/46`) already re-hashes the checkpoints against the end records, so a restored archive is checked automatically.

## Recommended approach: SHA-256 manifest in Git, lossless archive as a GitHub Release asset

1. **Manifest (committed).** `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` lists, for every file under `reports/phase9/corrected_rerun/`, `reports/phase9/corrected_rerun_report/` and `reports/phase9/dc_interval/`, its SHA-256, size and path, in `shasum -a 256` format, so `shasum -a 256 -c` checks it directly. It is generated once the campaign and its verification are final:
   ```bash
   cd /path/to/coad_transportability
   find reports/phase9/corrected_rerun reports/phase9/corrected_rerun_report reports/phase9/dc_interval -type f ! -name '*.lock' -print0 \
     | LC_ALL=C sort -z | xargs -0 shasum -a 256 > reports/phase9/ARCHIVE_MANIFEST_SHA256.txt
   ```
2. **Archive (not in Git).** A deterministic, uncompressed tar of exactly the listed files. Tar is lossless; gzip `-n` may be added because it is lossless too and omits timestamps:
   ```bash
   find reports/phase9/corrected_rerun reports/phase9/corrected_rerun_report reports/phase9/dc_interval -type f ! -name '*.lock' -print0 \
     | LC_ALL=C sort -z | tar --null -T - -cf phase9_campaign_artifacts.tar
   shasum -a 256 phase9_campaign_artifacts.tar > phase9_campaign_artifacts.tar.sha256
   ```
   The working copies are left in place; tar only reads them.
3. **Publish the archive** as an asset of a GitHub Release on the private repository (up to 2 GiB per asset; about 0.6 GB here). Release assets are stored outside Git history, so the repository stays small. For a public, citable copy, the same tar can go to Zenodo, which issues a DOI.
4. **Ignore only the oversized files in Git**, documented rather than silent. Add to `.gitignore`:
   ```
   # Phase 9 checkpoints (>100 MB at B_out = 750): archived as a release asset,
   # pinned by reports/phase9/ARCHIVE_MANIFEST_SHA256.txt
   reports/phase9/corrected_rerun/runs/*.jsonl
   reports/phase9/corrected_rerun/runs/*.lock
   ```
   All other campaign evidence (invocation records, gate evidence, logs, manifest, sizing file, verification reports) is small and is committed normally.
5. **Restore and check:**
   ```bash
   tar -xf phase9_campaign_artifacts.tar        # from the repository root
   shasum -a 256 -c reports/phase9/ARCHIVE_MANIFEST_SHA256.txt
   ```

## Archive handover

**Files that must stay outside ordinary Git tracking.** These five are the only files over 100 MB; every other campaign file is under 3 MB.

| File | Bytes | SHA-256 |
|---|---:|---|
| `reports/phase9/corrected_rerun/runs/run0.jsonl` | 112,951,341 | `0e27494ff875dcd7512355d213ee22f6bda2ce9a58595790b3a8b3b56b5dab15` |
| `reports/phase9/corrected_rerun/runs/run1.jsonl` | 112,957,953 | `9dfe640b32793b5e21d40fa4723c65f8455b15282f0ad443f7a26f645a9b7f43` |
| `reports/phase9/corrected_rerun/runs/run2.jsonl` | 112,960,361 | `0204da98f37ce4af5b5c7cdc61c114c602178f190bb25ca82dcba47b360e91a3` |
| `reports/phase9/corrected_rerun/runs/run3.jsonl` | 112,949,339 | `05274f0e5309454f9340e8a427304cd0d94a1c3538059784d85ba539fffc1e70` |
| `reports/phase9/corrected_rerun/runs/run4.jsonl` | 112,959,127 | `af7530077f5ce68ea874034d796a129327056afe20951451699d30559e58e326` |

The independent audit (`reports/phase9/independent_audit/`) confirmed for every file that:
- the hash equals the run's final end record;
- the 500-record prefix equals the 500-level end record;
- all 750 per-record hashes are valid.

**Manifest coverage** (52 files, about 590 MB; directory structure to preserve on restore, relative to the repository root):
```
reports/phase9/corrected_rerun/
    gate/          10 files   pre-execution gate evidence
    halts/         20 files   invocation start/end records (5 runs x 2 invocations)
    logs/           1 file
    manifest/       1 file    input_manifest.json (1ff24387...)
    runs/           5 files   the checkpoints above
    sizing/         1 file    level00500.json
    verification/   4 files   L00500 and L00750 reports
reports/phase9/corrected_rerun_report/   2 files   final interval results and tables
reports/phase9/dc_interval/              8 files   historical as-run checkpoints (14 MB; V1/V4 reference data)
```
`reports/phase9/independent_audit/` was created after the manifest. It is small and is committed normally, not archived.

**Recommended destination:** a GitHub Release asset on the private repository (`phase9_campaign_artifacts.tar`, about 0.6 GB, uncompressed and deterministic), with an optional Zenodo deposit for a citable DOI. **Not created; owner approval required.**

**`.gitignore` policy (proposed, not applied):** ignore only `reports/phase9/corrected_rerun/runs/*.jsonl` and `*.lock`, with the comment shown in step 4 above. Everything else, including the manifest, is tracked.

**Restore and verify:**
```bash
cd coad_transportability                       # repository root
tar -xf /path/to/phase9_campaign_artifacts.tar # recreates the directory structure above
shasum -a 256 -c reports/phase9/ARCHIVE_MANIFEST_SHA256.txt   # expect 52 lines ending in ": OK"
python3 -I -B scripts/48_phase9_independent_audit.py /tmp/p9_audit_check   # recomputes every reported number
```

## Alternatives considered

| Option | Verdict |
|---|---|
| Plain Git | Impossible: per-file limit of 100 MB |
| Git LFS | Works (2 GB per file), but GitHub's free LFS quota is 1 GB storage and 1 GB monthly bandwidth; each clone or CI run downloads about 0.6 GB. Acceptable for a private repository with few clones. |
| Splitting files into chunks | Rejected: the stored bytes would no longer be the files the pipeline hashed |
| Lossy compression or dropping fields | Rejected: destroys evidence |

## Status

- **Manifest generated (8 November 2026), after the PASS verification at B_out = 750:** `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt`, 52 files, about 590 MB in total. Its SHA-256 is `c45bb806…6b6f`, and `shasum -a 256 -c` reports all 52 OK. The five checkpoint hashes in it equal the `checkpoint_sha256` values in the runs' final end records.
- **Not done (owner decision):** the tar has not been built and nothing has been published or pushed. The `.gitignore` entry has not been added.
- **Owner decision required:** approve the Release-asset approach (or choose Git LFS) and the `.gitignore` entry above.
