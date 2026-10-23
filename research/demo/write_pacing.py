"""Locate the long live code-entry span for a short, explicit editing time lapse."""
import argparse
import json
from pathlib import Path
import pyte
from terminal import TerminalScreen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("recording", type=Path)
    args = parser.parse_args()
    evidence = json.loads((args.recording / "evidence.json").read_text())
    after = next(scene["start"] for scene in evidence["scenes"] if scene["scene"] == "consolidate the discount tests")
    screen = TerminalScreen(110, 32)
    stream = pyte.Stream(screen)
    start = end = None
    for line in (args.recording / "fast-tcp-demo.cast").read_text().splitlines()[1:]:
        seconds, kind, text = json.loads(line)
        if kind != "o":
            continue
        stream.feed(text)
        if seconds < after:
            continue
        if start is None and screen.display[0].lstrip().startswith("1 i"):
            start = seconds
        if start and any("assert discounted_price(price, percent) == expected" in row for row in screen.display):
            end = seconds
            break
    if start is None or end is None or not 15 < end - start < 50:
        raise SystemExit("Could not locate the live edit safely; inspect the capture before changing its pacing")
    result = {"intervals": [{"start": start, "end": end, "speed": 2.0,
                             "reason": "Live entry of the parameterized replacement tests; code-entry time lapse."}]}
    (args.recording / "pacing.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result))


if __name__ == "__main__":
    main()
