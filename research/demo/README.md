# Three-minute demonstration

[Watch the recording](recording/fast-tcp-demo.mp4) — **3:00**, 1920×1080, 20 fps, silent.

The video keeps the plain terminal/Vim presentation. A single 26-second explanation uses text, rectangles, and lines in the same palette to show character shingles, MinHash minima, candidate comparison, and the merged signature of selected tests. There are no title cards, sidebars, narration, or clock overlays.

| Video time | Action |
| --- | --- |
| 0:00–0:46 | Inspect the pricing project, initialize pytest, compare cold/warm partitions, inspect snapshot files |
| 0:46–1:12 | Five-character shingles → computed MinHash signatures → FAST-pw candidate ordering |
| 1:12–1:30 | Add two discount regression tests in Vim |
| 1:30–1:40 | Compare normal execution (50 passes before failure) with FAST (the regression fails first) |
| 1:40–1:49 | Fix the percentage divisor and run all 52 tests |
| 1:49–2:02 | Edit an existing tax test; its module's ten tests run first |
| 2:02–2:44 | Consolidate discount cases, verify the replacement tests, remove the old files, inspect the new/unchanged/deleted partition |
| 2:44–3:00 | Run FAST-one and finish with the unchanged 52-test suite |

## Evidence and scope

The terminal output comes from real subprocesses. `demo_observer.py` records collection order, setup/call/teardown outcomes, exit status, and source hashes without changing ordering. `recording/evidence.json` and the individual `.audit.json` files preserve the results, including the two intentionally stopped failing runs. The complete passing runs execute every collected test.

The pricing project is a constructed demonstration. A visible `conftest.py` fixture adds the same 30 ms delay per test used by the earlier short recording, so test progress is visible. The example's timings are not benchmark evidence.

The illustration uses three collected pytest node IDs and actual five-character shingle sets. `make_visual_data.py` computes the ten-component signatures with the package's `generate_signature`, checks their componentwise minima, and calls its `select_next_tests` for the candidate decisions. The displayed similarity estimates are 6/10 versus 2/10 relative to A, then 4/10 for B relative to the merged A+C signature. A is an explicit example seed. This illustrates candidate ranking; LSH candidate generation and the full suite's random initial choice are outside this small illustration. Hash values are shortened to hexadecimal prefixes for legibility; comparisons use full values.

The raw PTY capture is 170.775 seconds. One 28.924-second code-entry passage plays at 2×; the other terminal actions run approximately 1.015×. Together with the 26-second illustration, the result is exactly 180 seconds. No idle footage is inserted. The original timestamps and pacing decisions are preserved in the asciicast, `pacing.json`, and `render-manifest.json`.

## Reproduce

Requires Python 3.11, Git, Vim, and ffmpeg. From the repository root:

```bash
python3.11 -m venv output/video/demo-env
source output/video/demo-env/bin/activate
python -m pip install -e . pytest pexpect pillow pyte
export FAST_DEMO_VENV="$VIRTUAL_ENV"
export FAST_DEMO_PROJECT="$(mktemp -d)"
export FAST_DEMO_OUTPUT="$PWD/output/video/demo"
python research/demo/record_demo.py
python research/demo/make_visual_data.py "$FAST_DEMO_OUTPUT"
python research/demo/write_pacing.py "$FAST_DEMO_OUTPUT"
python research/demo/render_demo.py
```

Use a fresh project directory for each capture. `FAST_DEMO_PACE=0.01` runs a quick rehearsal with the same commands and assertions. A rehearsal is too short for the final renderer, which refuses to stretch a short capture to fill the duration.

The default font is macOS Menlo, matching the earlier recording. On another platform, set `FAST_DEMO_FONT` to a monospace font file before rendering. New capture runtimes may vary slightly; the renderer keeps the final duration at three minutes and records its pacing.
