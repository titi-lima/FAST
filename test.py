import time

import fast

if __name__ == "__main__":
    old_test_suite = [
        ("T001", "check login works"),
        ("T002", "check logout works"),  # del
        ("T003", "check password reset"),
        ("T004", "check user profile update"),
    ]
    new_test_suite = [
        ("T001", "check login works"),  # old
        ("T003", "check password reset"),  # old
        ("T004", "check user profile update"),  # old
        ("T005", "check email verification"),  # new
        ("T006", "check two-factor auth"),  # new
    ]

    start = time.perf_counter()
    fast.preparation(new_test_suite, old_test_suite)
    end = time.perf_counter() - start
    print(f"Preparation time: {end:.4f} seconds")

    start = time.perf_counter()
    prioritized_tests = fast.prioritization(new_test_suite, old_test_suite)
    end = time.perf_counter() - start
    print(f"Prioritization time: {end:.4f} seconds")
    print(prioritized_tests)

