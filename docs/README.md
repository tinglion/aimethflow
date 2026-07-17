# AIMethFlow Documentation

Detailed documentation for the AIMethFlow plasma cfDNA methylation pipeline.
For a quick overview, start with the [project README](../README.md).

## Getting started

| Document | Purpose |
|---|---|
| [installation.md](installation.md) | Conda environment setup, reference genome indexing, auxiliary data files |
| [../DOCKER.md](../DOCKER.md) | Docker-based deployment (recommended for production) |
| [quickstart.md](quickstart.md) | Web UI walkthrough, full-pipeline command, step-by-step examples, shell-script template |

## Reference

| Document | Purpose |
|---|---|
| [pipeline.md](pipeline.md) | 20-step pipeline reference (plasma and tissue modes); detailed Step 19 cross-validation section |
| [configuration.md](configuration.md) | `config.toml` field-by-field reference |
| [outputs.md](outputs.md) | Output file catalogue and sample Excel format |
| [resources.md](resources.md) | Hardware sizing, peak memory, disk-usage measurements |

## Methodology and operations

| Document | Purpose |
|---|---|
| [methodology.md](methodology.md) | Library-preparation prerequisites, methodological improvements over prior approaches, UMI-aware deduplication deep-dive with quantitative audits |
| [tools.md](tools.md) | Helper scripts: FastQ inspection, panel BED derivation, UMI-rescue audits, cleanup commands |

## Related project memos

Internal design notes and assay-improvement memos live in this same directory; they are not part of the user-facing documentation tree:

- `sting_血浆分析参考组织分析改进建议.md` — design memo on tissue-vs-plasma reference choices
