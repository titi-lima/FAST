"""Compute the animation's shingles, minima, and FAST-pw candidate decisions."""
import argparse
import json
from pathlib import Path
import numpy as np
from fast_tcp import fast


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("recording", type=Path)
    args = parser.parse_args()
    collected = json.loads((args.recording / "cold.audit.json").read_text())["collected"]
    names = {
        "A": next(node for node in collected if node.endswith("test_shipping[49-False-5]")),
        "B": next(node for node in collected if node.endswith("test_shipping[50-False-0]")),
        "C": next(node for node in collected if node.endswith("test_tax[100-20-20]")),
    }
    fast.DEFAULT_K = 5
    fast.DEFAULT_R = 1
    fast.DEFAULT_B = 10
    fast.DEFAULT_H = 10
    fast.DEFAULT_ALG = "FAST-pw"
    signatures = {key: fast.generate_signature(name) for key, name in names.items()}
    # Visualize the real candidate-ranking function. LSH candidate generation is
    # outside this small illustration; A is an explicit example seed, not a
    # claim that this test is the random first choice for the whole suite.
    candidates = ["B", "C"]
    cumulative = signatures["A"].copy()
    order = ["A"]
    decisions = []
    while candidates:
        scores = {key: cumulative.jaccard(signatures[key]) for key in candidates}
        chosen = fast.select_next_tests(candidates, signatures, cumulative)[0]
        decisions.append({"scores": scores, "selected": chosen})
        order.append(chosen)
        cumulative.merge(signatures[chosen])
        candidates.remove(chosen)
    assert order == ["A", "C", "B"], order
    shingles = sorted(fast.k_shingles(names["A"]))
    per_shingle = {shingle: fast.generate_signature(shingle).hashvalues.tolist() for shingle in shingles}
    # A five-character input contains exactly one five-character shingle.
    assert np.array_equal(np.min(list(per_shingle.values()), axis=0), signatures["A"].hashvalues)
    winner = min(shingles, key=lambda shingle: per_shingle[shingle][0])
    shipping_start = names["A"].index("shipping")
    windows = [names["A"][offset:offset + 5] for offset in range(shipping_start, shipping_start + 4)]
    examples = [winner] + [shingle for shingle in windows if shingle != winner][:3]
    merged_ac = signatures["A"].copy()
    merged_ac.merge(signatures["C"])
    data = {
        "k": 5, "permutations": 10,
        "scope": "FAST-pw comparison of supplied candidates; LSH filtering omitted; A is an explicit example seed.",
        "tests": {key: {"nodeid": names[key], "display": names[key].split("::")[1],
                         "hashvalues": [int(value) for value in signatures[key].hashvalues]}
                  for key in names},
        "windows": windows, "window_start": shipping_start,
        "shingle_examples": [{"text": value, "h1": per_shingle[value][0]} for value in examples],
        "minimum_shingle": winner, "decisions": decisions, "order": order,
        "cumulative_ac": [int(value) for value in merged_ac.hashvalues],
        "verified": True,
    }
    (args.recording / "visual-data.json").write_text(json.dumps(data, indent=2))
    print(json.dumps({"candidate_decisions": decisions, "order": order, "verified": True}), flush=True)


if __name__ == "__main__":
    main()
