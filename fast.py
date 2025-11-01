import json
import math
import os
import pickle
import random

import xxhash
from datasketch import MinHash, MinHashLSH

###############################################################################
# CONFIGURATION PARAMETERS


# Load configuration from JSON file
CONFIG_PATH = "fast_config.json"
if os.path.exists(CONFIG_PATH):
    try:
        with open(CONFIG_PATH) as f:
            config = json.load(f)
    except (json.JSONDecodeError, OSError):
        config = {}
else:
    config = {}

# Default configuration parameters. Change in the JSON file only!
# If the file is missing or malformed, defaults are used
SIGNATURE_FOLDER = config.get("signature_folder", "signatures")  # Folder of MinHash signatures
DEFAULT_K = config.get("k", 5)  # k-shingles parameter
DEFAULT_R = config.get("r", 1)  # lsh: number of rows
DEFAULT_B = config.get("b", 10)  # lsh: number of bands
DEFAULT_H = DEFAULT_R * DEFAULT_B  # minhash: number of permutations = r*b
DEFAULT_ALG = config.get("alg", "FAST-pw")  # FAST prioritization algorithm
DEFAULT_BUDGET = config.get("budget", 0)  # prioritization budget (0 = all tests)

assert SIGNATURE_FOLDER, "signature_folder must be non-empty"
assert DEFAULT_K > 0, "k must be positive"
assert DEFAULT_R > 0, "r must be positive"
assert DEFAULT_B > 0, "b must be positive"
assert DEFAULT_H == DEFAULT_R * DEFAULT_B, "h must equal r * b"
assert DEFAULT_ALG in {"FAST-pw", "FAST-one", "FAST-log", "FAST-sqrt", "FAST-all"}, "unknown alg"
assert DEFAULT_BUDGET >= 0, "budget must be non-negative"


###############################################################################
# SIGNATURES


def k_shingles(document):
    """
    Return the set of k-shingles (contiguous substrings of length k) for `text`.

    Parameters
    - text: Input string from which to extract shingles.

    Returns
    - A set of unique k-length substrings from `text`.
    """
    k = DEFAULT_K

    return {document[i : i + k] for i in range(len(document) - k + 1)}


def generate_signature(document):
    """
    Generate a datasketch.MinHash signature for `document`.

    Parameters
    - document: Input text to be fingerprinted.

    Returns
    - A datasketch.MinHash instance representing the document's MinHash signature.

    Implementation details
    - Uses xxhash.xxh64_intdigest as the MinHash hashfunc for consistent 64-bit hashing.
    - Each shingle is encoded as UTF-8 bytes before updating the MinHash.
    """
    h = DEFAULT_H

    shingles = k_shingles(document)
    m = MinHash(num_perm=h, hashfunc=xxhash.xxh64_intdigest)
    for s in shingles:
        m.update(s.encode("utf8"))
    return m


def write_signature(obj, path):
    """
    Serialize `obj` to `path` using pickle in binary mode.

    Parameters
    - obj: Any picklable Python object (e.g., a MinHash instance).
    - path: Filesystem path where the pickle will be written.
    """
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def read_signature(path):
    """
    Deserialize and return a pickled object from `path`.

    Parameters
    - path: Filesystem path to a pickle file.

    Returns
    - The deserialized Python object.
    """
    with open(path, "rb") as f:
        return pickle.load(f)


def load_signatures(new_test_suite):
    """
    Load persisted (pickle serialization) MinHash signatures for tests in the new test suite.

    Parameters:
    - new_test_suite: iterable of (t_id, test_content).

    Returns:
    - signatures: dict mapping t_id -> MinHash signature object.
    """
    path = SIGNATURE_FOLDER

    signatures = {}
    for t_id, _ in new_test_suite:
        sig = read_signature(os.path.join(path, f"{t_id}.pkl"))
        signatures[t_id] = sig
    return signatures


###############################################################################
# PARTITION TEST SUITE


def partition_test_suite(old_test_suite, new_test_suite):
    """
    Partition tests into unchanged, deleted, and new sets by comparing
    (test_id, hash(test)) pairs from old and new test suites.

    Parameters:
    - old_test_suite: iterable of (t_id, test_content) representing the previous suite.
    - new_test_suite: iterable of (t_id, test_content) representing the current suite.

    Returns:
    - old_tests: set of t_id present in both suites with identical content.
    - del_tests: set of t_id present in old suite but not in new suite (or changed).
    - new_tests: set of t_id present in new suite but not in old suite (or changed).
    """
    old_hash = {(t_id, xxhash.xxh64_intdigest(t.encode("utf8"))) for t_id, t in old_test_suite}
    new_hash = {(t_id, xxhash.xxh64_intdigest(t.encode("utf8"))) for t_id, t in new_test_suite}

    new_tests = {t_id for t_id, _ in new_hash - old_hash}
    old_tests = {t_id for t_id, _ in old_hash & new_hash}
    del_tests = {t_id for t_id, _ in old_hash - new_hash}

    return new_tests, old_tests, del_tests


###############################################################################
# FAST PREPARATION


def preparation(new_test_suite, old_test_suite):
    """
    Ensure the signature storage directory exists and synchronize persisted
    signatures with the new test suite.

    Actions:
    - Remove persisted signatures for tests deleted from the suite.
    - Generate and store signatures for tests that do not yet have persisted files.

    Parameters:
    - new_test_suite: iterable of (t_id, test_content).
    - old_test_suite: iterable of (t_id, test_content) from the previous run.
    """
    path = SIGNATURE_FOLDER

    _, _, del_tests = partition_test_suite(old_test_suite, new_test_suite)

    os.makedirs(path, exist_ok=True)
    for t_id in del_tests:
        if os.path.exists(os.path.join(path, f"{t_id}.pkl")):
            os.remove(os.path.join(path, f"{t_id}.pkl"))
    for t_id, t in new_test_suite:
        if not os.path.exists(os.path.join(path, f"{t_id}.pkl")):
            sig = generate_signature(t)
            write_signature(sig, os.path.join(path, f"{t_id}.pkl"))


###############################################################################
# FAST PRIORITIZATION


def lsh_buckets(test_suite, remaining_tests, signatures):
    """
    Build an LSH index (MinHashLSH) for the tests currently remaining.

    Parameters:
    - test_suite: set of t_id.
    - remaining_tests: set of t_id to include in the index.
    - signatures: dict mapping t_id -> MinHash signature.

    Returns:
    - lsh: MinHashLSH instance with remaining tests inserted.
    """
    r = DEFAULT_R
    b = DEFAULT_B
    h = DEFAULT_H

    lsh = MinHashLSH(num_perm=h, params=(b, r))
    for t_id in test_suite:
        sig = signatures[t_id]
        if t_id in remaining_tests:
            lsh.insert(t_id, sig)
    return lsh


def cumulative_signature(prioritized_test_suite, signatures):
    """
    Build a cumulative MinHash signature that merges signatures of tests
    that have already been prioritized (i.e., not in remaining_tests).

    Parameters:
    - test_suite: set of t_id.
    - remaining_tests: set of t_id still to be prioritized.
    - signatures: dict mapping t_id -> MinHash signature.

    Returns:
    - cumulative_sig: MinHash instance representing the merge of selected tests.
    """
    h = DEFAULT_H

    cumulative_sig = MinHash(num_perm=h, hashfunc=xxhash.xxh64_intdigest)
    for t_id in prioritized_test_suite:
        sig = signatures[t_id]
        cumulative_sig.merge(sig)
    return cumulative_sig


def generate_candidates(lsh, remaining_tests, cumulative_sig):
    """
    Generate a candidate set of tests to select next.

    The function queries the LSH with the cumulative signature to find
    tests similar to the already selected set, then returns remaining_tests
    minus those similar candidates. If this yields an empty candidate set,
    it clears the cumulative signature and retries once; if still empty,
    all remaining tests are returned.

    Parameters:
    - lsh: MinHashLSH index over remaining tests.
    - remaining_tests: set of t_id still to be prioritized.
    - cumulative_sig: MinHash representing merged selected signatures.

    Returns:
    - list of candidate t_id.
    """
    similar_candidates = set(lsh.query(cumulative_sig))
    candidates = remaining_tests - similar_candidates
    if not candidates:
        cumulative_sig.clear()
        similar_candidates = set(lsh.query(cumulative_sig))
        candidates = remaining_tests - similar_candidates
        if not candidates:
            candidates = remaining_tests
    return list(candidates)


def candidate_set_size(x):
    """
    Compute the size of the candidate subset to pick based on the algorithm.

    Parameters:
    - x: integer current number of candidates.

    Returns:
    - integer in {1,...,x} representing how many candidates to select.

    Raises:
    - ValueError if alg is not present in the `algs` mapping.
    """
    alg = DEFAULT_ALG

    # Mapping of prioritization algorithms to their candidate set size functions
    algs = {
        "FAST-one": lambda x: 1,
        "FAST-log": lambda x: math.log2(x),
        "FAST-sqrt": lambda x: math.sqrt(x),
        "FAST-all": lambda x: x,
    }
    try:
        set_size = int(algs[alg](x))
    except KeyError:
        raise ValueError(f"Unknown algorithm: {alg}")
    return max(1, min(set_size, x))


def select_next_tests(candidates, signatures, cumulative_sig):
    """
    Select the next test(s) to add to the prioritized suite.

    Behavior:
    - If alg == "FAST-pw": pick the single test in `candidates` that is
      most dissimilar to the cumulative signature (minimizes Jaccard similarity).
    - Otherwise: treat `alg` as one of the size-determining functions and
      randomly sample that many tests from `candidates`, then shuffle.

    Parameters:
    - candidates: list of candidate t_id.
    - signatures: dict mapping t_id -> MinHash signature.
    - cumulative_sig: MinHash representing merged selected signatures.

    Returns:
    - next_tests: list of t_id selected next (shuffled).
    """
    alg = DEFAULT_ALG

    next_tests = []
    if alg == "FAST-pw":
        most_dissimilar_test = min(
            candidates, key=lambda t_id: cumulative_sig.jaccard(signatures[t_id])
        )
        next_tests.append(most_dissimilar_test)
    else:  # fast-f
        sel_size = candidate_set_size(len(candidates))
        next_tests = random.sample(candidates, sel_size)
    random.shuffle(next_tests)
    return next_tests


def fast_alg(test_suite, signatures, prioritized_test_suite, budget):
    """
    Produce an ordering for tests in `test_suite` using the FAST algorithm.

    This function orders only the tests present in `test_suite`. The cumulative
    signature used to define dissimilarity is seeded with `prioritized_test_suite`
    (tests already selected before this call). The returned list contains only
    tests from `test_suite` (it does not include the seed items).

    Behavior summary:
    - If no seed is provided, one random test from `test_suite` is chosen first.
    - An LSH index is built over remaining tests' MinHash signatures.
    - Repeatedly generate candidate sets (tests dissimilar to the cumulative
      signature), select next tests via `select_next_tests`, merge their
      signatures into the cumulative signature, and append them to the local
      ordering until `budget` items have been selected or the set is exhausted.
    - The LSH is rebuilt when the remaining set shrinks below a halving threshold.

    Parameters:
    - test_suite: iterable or set of t_id to order.
    - signatures: dict mapping t_id -> MinHash signature (loaded via load_signatures).
    - prioritized_test_suite: sequence of already-selected t_id used only to seed
      the cumulative signature.

    Returns:
    - list of t_id from `test_suite` ordered according to FAST (length <= budget).
    """
    prioritized_local = []
    remaining_tests = set(test_suite)
    if not prioritized_test_suite:
        first_test = random.choice(list(remaining_tests))
        prioritized_local.append(first_test)
        remaining_tests.remove(first_test)

    lsh = lsh_buckets(test_suite, remaining_tests, signatures)
    prioritized_global = prioritized_test_suite + prioritized_local
    cumulative_sig = cumulative_signature(prioritized_global, signatures)

    threshold = len(test_suite)
    while remaining_tests and len(prioritized_local) < budget:
        # recompute lsh bucket every time we halve the number of remaining tests
        # this makes the computation of candidates more efficient
        if len(remaining_tests) < threshold:
            lsh = lsh_buckets(test_suite, remaining_tests, signatures)
            threshold /= 2

        candidates = generate_candidates(lsh, remaining_tests, cumulative_sig)
        next_tests = select_next_tests(candidates, signatures, cumulative_sig)
        for t_id in next_tests:
            cumulative_sig.merge(signatures[t_id])
            prioritized_local.append(t_id)
            remaining_tests.remove(t_id)

    return prioritized_local[:budget]


def prioritization(new_test_suite, old_test_suite):
    """
    High-level FAST prioritization entry point.

    Algorithm:
    1. Partition the suite into new and unchanged/old tests.
    2. Load stored MinHash signatures for the current suite.
    3. Order new tests first by calling `fast_alg` (so new tests get higher priority).
    4. If more tests are needed to meet the budget, order old (unchanged) tests next,
       seeding the second phase with signatures of already-selected new tests so
       duplicates and overly-similar tests are avoided.
    5. Return the prioritized list truncated to `budget`.

    Parameters:
    - new_test_suite: iterable of (t_id, test_content) for the current run.
    - old_test_suite: iterable of (t_id, test_content) from the previous run.

    Returns:
    - prioritized_test_suite: list of t_id ordered by priority, length <= budget.
    """
    budget = DEFAULT_BUDGET

    if budget == 0:
        budget = len(new_test_suite)

    signatures = load_signatures(new_test_suite)
    new_tests, old_tests, _ = partition_test_suite(old_test_suite, new_test_suite)
    prioritized_test_suite = []
    if new_tests:
        budget_new = budget
        prioritized_new = fast_alg(new_tests, signatures, prioritized_test_suite, budget_new)
        prioritized_test_suite.extend(prioritized_new)
    if len(prioritized_test_suite) < budget:
        budget_old = budget - len(prioritized_test_suite)
        prioritized_old = fast_alg(old_tests, signatures, prioritized_test_suite, budget_old)
        prioritized_test_suite.extend(prioritized_old)
    return prioritized_test_suite[:budget]


