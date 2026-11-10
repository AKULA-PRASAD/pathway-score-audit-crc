# Commit ID map (history rewrite, 9 November 2026)

On 9 November 2026 the five Phase 10 commits were rebuilt to remove tool co-author lines from their commit messages. **Files, authors and dates are unchanged**: each rebuilt commit has exactly the same file tree as the original. Everything up to and including `v1.0.0` (`9cba586`) keeps its original IDs, as do the `phase9-corrected-campaign-2026-11-08` tag and the `extensions-tool-portfolio` branch.

| Original | Rebuilt | Commit |
|---|---|---|
| `125667a26a173925e72a4224639ac498ddebc7b4` | `007ff3b` | Phase 10: register and freeze Track B before TCGA unblinding |
| `2c4ef3aaefa24af997d8678c0aa8fb306357090d` | `97bdccf` | Phase 10: re-freeze after evaluation-metric fix (declared deviation 1) |
| `0877d2c5ca4155085ab5d4ddc78c82ab847f8904` | `9935802` | Phase 10: record non-informative HiSeqV2 evaluation; register GA evaluation (deviation 2) |
| `c7e7856` | `636c05a` | Phase 10: results, documentation, website ML section and resume text |
| `2ccd188` | `04e10cf` | Merge Phase 10: leakage-controlled ML upgrade (PR #5) |

## Records that cite the original IDs

These sealed records are left unchanged, because editing them would break their own integrity checks. Read their commit IDs through the table above.

| File | Field | Original → rebuilt |
|---|---|---|
| `reports/phase10/tcga/UNBLINDING_APPROVAL.json` | `freeze_commit` | `2c4ef3a` → `97bdccf` |
| `reports/phase10/tcga_ga/UNBLINDING_APPROVAL.json` | `registration_commit` | `0877d2c` → `9935802` |
| `reports/phase10/tcga/freeze_manifest.json` | `git_head` | `125667a` → `007ff3b` |
| `reports/phase10/mmr/summary.json` | `git_head` | `0877d2c` → `9935802` |

The order of events is unchanged: the freeze and the GA registration were committed before the TCGA labels were read. The original commits remain viewable on pull request #5.
