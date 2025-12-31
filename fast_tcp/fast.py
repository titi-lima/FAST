import json
import math
import os
import pickle
import random
from typing import Dict, Iterable, List, Set, Tuple, Union

import xxhash
from datasketch import MinHash, MinHashLSH  # type: ignore

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
SIGNATURE_FOLDER = config.get(
    "signature_folder", "signatures"
)  # Folder of MinHash signatures
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
assert DEFAULT_ALG in {
    "FAST-pw",
    "FAST-one",
    "FAST-log",
    "FAST-sqrt",
    "FAST-all",
}, "unknown alg"
assert DEFAULT_BUDGET >= 0, "budget must be non-negative"


###############################################################################
# SIGNATURES

# In-memory cache for signatures to avoid repeated disk I/O
# Format: {test_id: (content_hash, MinHash)}
_SIGNATURE_CACHE: Dict[str, Tuple[int, MinHash]] = {}
_SIGNATURE_CACHE_PATH: str = ""


def _content_hash(content: str) -> int:
    """Compute a fast hash of test content for change detection."""
    return xxhash.xxh64_intdigest(content.encode("utf-8"))


def k_shingles(document: str) -> Set[str]:
    """
    Return the set of k-shingles (contiguous substrings of length k) for `text`.

    Parameters
    - text: Input string from which to extract shingles.

    Returns
    - A set of unique k-length substrings from `text`.
    """
    k = DEFAULT_K

    return {document[i : i + k] for i in range(len(document) - k + 1)}


def generate_signature(document: str) -> MinHash:
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


def _get_consolidated_signature_path() -> str:
    """Get the path to the consolidated signature file."""
    return os.path.join(SIGNATURE_FOLDER, "signatures.pkl")


def _load_signature_cache() -> Dict[str, Tuple[int, MinHash]]:
    """Load the consolidated signature cache from disk.
    
    Returns dict mapping test_id -> (content_hash, MinHash).
    """
    global _SIGNATURE_CACHE, _SIGNATURE_CACHE_PATH
    
    consolidated_path = _get_consolidated_signature_path()
    
    # Return cached data if already loaded for this path
    if _SIGNATURE_CACHE_PATH == consolidated_path and _SIGNATURE_CACHE:
        return _SIGNATURE_CACHE
    
    _SIGNATURE_CACHE_PATH = consolidated_path
    
    if os.path.exists(consolidated_path):
        try:
            with open(consolidated_path, "rb") as f:
                _SIGNATURE_CACHE = pickle.load(f)
        except (pickle.PickleError, OSError, EOFError):
            _SIGNATURE_CACHE = {}
    else:
        _SIGNATURE_CACHE = {}
    
    return _SIGNATURE_CACHE


def _save_signature_cache(signatures: Dict[str, Tuple[int, MinHash]]) -> None:
    """Save the consolidated signature cache to disk."""
    global _SIGNATURE_CACHE, _SIGNATURE_CACHE_PATH
    
    consolidated_path = _get_consolidated_signature_path()
    os.makedirs(SIGNATURE_FOLDER, exist_ok=True)
    
    with open(consolidated_path, "wb") as f:
        pickle.dump(signatures, f)
    
    _SIGNATURE_CACHE = signatures
    _SIGNATURE_CACHE_PATH = consolidated_path


def load_signatures(new_test_suite: List[Tuple[str, str]]) -> Dict[str, MinHash]:
    """
    Load persisted MinHash signatures for tests in the new test suite.

    Uses a consolidated pickle file for all signatures instead of individual
    files per test, reducing I/O overhead from O(n) to O(1).

    Parameters:
    - new_test_suite: iterable of (t_id, test_content).

    Returns:
    - signatures: dict mapping t_id -> MinHash signature object.
    """
    cache = _load_signature_cache()

    signatures = {}
    for t_id, _ in new_test_suite:
        if t_id in cache:
            # Extract just the MinHash from (content_hash, MinHash) tuple
            signatures[t_id] = cache[t_id][1]

    return signatures


###############################################################################
# FAST PREPARATION


def preparation(test_suite: List[Tuple[str, str]], del_tests: Set[str]) -> None:
    """
    Ensure the signature storage directory exists and synchronize persisted
    signatures with the new test suite.

    Actions:
    - Remove persisted signatures for tests deleted from the suite.
    - Generate and store signatures for new tests or tests whose content changed.

    Parameters:
    - test_suite: iterable of (t_id, test_content).
    - del_tests: set of t_id representing tests deleted from the previous suite.

    Uses a single consolidated pickle file for all signatures instead of
    individual files per test, reducing I/O overhead from O(n) to O(1).
    """
    # Load existing signature cache
    cache = _load_signature_cache()
    modified = False

    # Remove signatures for deleted tests
    for t_id in del_tests:
        if t_id in cache:
            del cache[t_id]
            modified = True

    # Generate signatures for new tests or tests whose content changed
    for t_id, t in test_suite:
        current_hash = _content_hash(t)
        if t_id not in cache:
            # New test - generate signature
            cache[t_id] = (current_hash, generate_signature(t))
            modified = True
        else:
            # Existing test - check if content changed
            cached_hash, _ = cache[t_id]
            if cached_hash != current_hash:
                # Content changed - regenerate signature
                cache[t_id] = (current_hash, generate_signature(t))
                modified = True

    # Save consolidated cache if anything changed
    if modified:
        _save_signature_cache(cache)


###############################################################################
# FAST PRIORITIZATION


def lsh_buckets(
    test_suite: Iterable[str],
    remaining_tests: Set[str],
    signatures: Dict[str, MinHash],
) -> MinHashLSH:
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


def cumulative_signature(
    prioritized_test_suite: Iterable[str], signatures: Dict[str, MinHash]
) -> MinHash:
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


def generate_candidates(
    lsh: MinHashLSH, remaining_tests: Set[str], cumulative_sig: MinHash
) -> List[str]:
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


def candidate_set_size(x: int) -> int:
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


def select_next_tests(
    candidates: List[str], signatures: Dict[str, MinHash], cumulative_sig: MinHash
) -> List[str]:
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


def fast_alg(
    test_suite: Set[str],
    signatures: Dict[str, MinHash],
    prioritized_test_suite: List[str],
    budget: int,
) -> List[str]:
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

    threshold: Union[int, float] = len(test_suite)
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


def prioritization(
    test_suite: List[Tuple[str, str]],
    new_tests: Set[str],
    old_tests: Set[str],
) -> List[str]:
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
    - test_suite: iterable of (t_id, test_content) for the current run.
    - new_tests: set of t_id representing tests new to the suite.
    - old_tests: set of t_id representing unchanged tests in the suite.

    Returns:
    - prioritized_test_suite: list of t_id ordered by priority, length <= budget.
    """
    budget = DEFAULT_BUDGET

    if budget == 0:
        budget = len(test_suite)

    signatures = load_signatures(test_suite)
    prioritized_test_suite: List[str] = []

    if new_tests:
        budget_new = budget
        prioritized_new = fast_alg(
            new_tests, signatures, prioritized_test_suite, budget_new
        )
        prioritized_test_suite.extend(prioritized_new)

    if len(prioritized_test_suite) < budget:
        budget_old = budget - len(prioritized_test_suite)
        prioritized_old = fast_alg(
            old_tests, signatures, prioritized_test_suite, budget_old
        )
        prioritized_test_suite.extend(prioritized_old)

    return prioritized_test_suite[:budget]
