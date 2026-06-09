# Docker Deployment

This is the recommended deployment path for production / shared servers. For local development with Conda see [docs/installation.md](docs/installation.md).

## Usage

**Prerequisites**: Docker ≥ 24, Docker Compose Plugin (`docker compose` command available)

```bash
# 1. Copy and fill in paths
cp .env.example .env
# Edit .env: set the host machine's data / reference genome / results directories

# 2. Build the image (first time ~10–20 minutes; cached on subsequent runs)
docker compose --progress=plain build

# 3. Start the persistent container
docker compose up -d
# Submit an analysis job
docker compose exec analysis \
    python /app/src/scripts/run_analysis.py --steps 1-19
# Enter the container shell
docker compose exec analysis bash

# One-shot run (removed after completion)
docker compose --profile run run --rm analysis-run \
    python /app/src/scripts/run_analysis.py --steps 1-19

# 4. Start JupyterLab to browse results step by step
docker compose --profile jupyter up jupyter
#    Open http://localhost:8888 in your browser
```

The Web UI is also auto-started in the persistent container and listens on port 9999 — see [docs/quickstart.md](docs/quickstart.md#web-ui-recommended).

## Protected core via mounted wheel (pre-publication only)

Steps 14 (feature selection), 15 (diagnostic classifiers), and 19 (external validation) depend on the `aimethflow-core` package, temporarily held private pending iMeta 2026 publication. See the *Status — pre-publication* section in [README.md](README.md) for context.

The container ships with stubs for these 3 modules. The pipeline imports cleanly and Steps 1-13, 16-18, 20 run end-to-end without any extra setup. Attempting Step 14 / 15 / 19 without the core installed raises a clean `ImportError`.

To enable the full pipeline, iMeta-assigned reviewers receive a single bytecode-only wheel `aimethflow_core-<version>-py310-none-any.whl` (~48 KB) from the corresponding author through the iMeta editorial system. Setup:

```bash
# 1. Place the wheel in a directory called `private/` next to docker-compose.yml
mkdir -p private
cp ~/Downloads/aimethflow_core-0.1.0-py310-none-any.whl private/

# 2. (Re-)start the container; the entrypoint auto-installs the wheel on boot
docker compose up -d
# Or for an already-running container:
docker compose restart analysis
```

On startup `docker-entrypoint.sh` checks for `/private/aimethflow_core-*.whl` (mounted from `./private/` on the host), pip-installs it if found, and exec's the original command. Subsequent restarts skip the install if the package is already present.

Without the wheel mounted, the container still starts and runs the public pipeline normally — the entrypoint logs a one-line notice that Steps 14, 15, 19 are unavailable.

The `private/` directory is gitignored. Do not commit the wheel.

## Design decisions

| Item | Decision |
|------|------|
| Base image | `continuumio/miniconda3` |
| Bioinformatics tools | conda (bioconda + conda-forge), consistent with the server environment |
| Python dependencies | conda env `meth` (Python 3.10.14) + pip |
| GPU | Not enabled |
| Reference genome | Mounted externally, not bundled |
| `data` / `temp` / `result` | Mounted externally, not bundled |
| Jupyter | Independent service in docker-compose, sharing volumes with the analysis container |

## Image layer structure

```text
continuumio/miniconda3
  ├─ [conda] mamba (solver acceleration)
  ├─ [mamba] create -n meth python=3.10.14    ~400 MB
  ├─ [mamba] bioconda bioinformatics tools     ~800 MB
  │          fastp / bwameth / samtools /
  │          bedtools / methyldackel / umi_tools / pysam
  ├─ [pip]   requirements.txt + jupyterlab    ~500 MB
  └─ [COPY]  src/                              ~10 MB
             (data / temp / result mounted at runtime)
```

## Tool installation method

### conda (bioconda + conda-forge)

| Tool | Version | Purpose |
|------|------|------|
| fastp | 0.23.4 | Quality control + trimming |
| bwameth | latest | Bisulfite alignment |
| samtools | latest | BAM processing |
| bedtools | latest | Region filtering |
| methyldackel | latest | Methylation extraction |
| umi_tools | latest | UMI-aware deduplication |
| pysam | 0.22.0 | Python BAM reader |

### pip (inside the conda env)

```text
jupyterlab ipykernel   # Jupyter environment
-r requirements.txt    # Project Python dependencies
```

## Mounted volumes

| Host path (example) | Container path | Read-only | Description |
| --- | --- | --- | --- |
| `./data` | `/app/data` | yes | panel BED, sample Excel, candidate regions |
| `./config.toml` | `/app/config.toml` | yes | Project configuration |
| `/data/meth/OC_tissue` | `/mnt/OC_tissue` | yes | Tissue bedGraphs |
| `/data/meth/OC_plasma` | `/mnt/OC_plasma` | yes | Plasma BAMs / bedGraphs |
| `/data/meth/ref` | `/mnt/ref` | yes | hg38 + bwa-meth index |
| `/data/meth/temp` | `/mnt/temp` | | Intermediate files |
| `/data/meth/result` | `/mnt/result` | | Analysis results |

Update paths in `config.toml` (or use `config.docker.toml`) to the container-internal paths, for example:

```toml
tissue_dir         = "/mnt/OC_tissue"
plasma_dir         = "/mnt/OC_plasma"
reference_genome   = "/mnt/ref/hg38.fa"

[directories]
temp_dir   = "/mnt/temp"
result_dir = "/mnt/result"
```

## Ports

| Port | Service |
|------|------|
| 8888 | JupyterLab |
| 9999 | Web UI (FastAPI) |

## docker-compose service plan

```yaml
services:
  analysis:          # Full pipeline (python run_analysis.py)
    build: .
    volumes: [data, config, OC_tissue, OC_plasma, ref, temp, result]
    command: python /app/src/scripts/run_analysis.py --steps 1-19

  jupyter:           # Step-by-step result inspection
    build: .
    ports: ["8888:8888"]
    volumes: [data, config, OC_tissue, OC_plasma, ref, temp, result, notebooks]
    command: jupyter lab --ip=0.0.0.0 --no-browser --allow-root
```

## Environment variables

| Variable | Default | Description |
|------|------|------|
| `NCBI_API_KEY` | *(empty)* | Raises the PubMed rate limit from 3 to 10 req/s (optional) |

## File inventory (to be created)

```text
meth/
├── Dockerfile            # Image build
├── docker-compose.yml    # One-command launch (analysis + jupyter)
├── .dockerignore         # Excludes data/, temp/, result/, __pycache__
├── DOCKER.md             # This file
└── notebooks/            # Jupyter step-by-step analysis notebooks
    ├── 01_qc.ipynb
    ├── ...
    └── 19_report.ipynb
```
