# Research artifacts

Start with [replication/README.md](replication/README.md) for the experiment commands, pinned project revisions and dependencies, recorded full-suite reports, and measured figure renderers. Everything needed for those experiments is in this repository; no Overleaf files are needed.

| Material | Status |
| --- | --- |
| [replication/](replication/README.md) | Verified real-project case study, executed synthetic workload, standalone snapshot experiment, recorded evidence, and reproduction instructions |
| [demo/](demo/README.md) | Three-minute terminal/Vim demonstration with a computed shingles/MinHash illustration, recording scripts, and execution evidence |
| [snapshot_cache_growth.ipynb](snapshot_cache_growth.ipynb) | Original snapshot experiment, including saved notebook output; standalone version in `replication/measure_snapshot.py` |
| [real_project_case_study.ipynb](real_project_case_study.ipynb) | Preliminary case study; superseded by the complete-suite measurements in `replication/data/` |
| [incremental_benefit_benchmark.ipynb](incremental_benefit_benchmark.ipynb) | Exploratory collection-only benchmark; does not measure complete test execution |
| [multiframework_benchmark.ipynb](multiframework_benchmark.ipynb) | Exploratory integration runs, including unsuccessful runs; not the verified case study |
| [profile_benchmark.py](profile_benchmark.py) | Separate profiling workload |
| [figures/generate_figures.py](figures/generate_figures.py) | Legacy plots with preliminary values and constructed timings; requires `--legacy`. Use the measured renderers in `replication/` for the verified results |
| [paper.tex](paper.tex), [poster.tex](poster.tex) | Local manuscript and poster sources; the manuscript may lag the version being edited in Overleaf |

The replication instructions describe the scope of each experiment. Runtime measurements establish execution cost and outcome preservation for the measured suites; they do not establish fault-detection effectiveness.
