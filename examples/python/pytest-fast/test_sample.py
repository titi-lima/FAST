import time


def test_c():
    time.sleep(0.01)
    assert 1 + 1 == 2


def test_a():
    time.sleep(0.01)
    assert "a".upper() == "A"


def test_b():
    time.sleep(0.01)
    assert sorted([3, 1, 2]) == [1, 2, 3]
