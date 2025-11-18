import time


def test_d():
    time.sleep(0.01)
    assert 1 + 1 == 2


def test_e():
    time.sleep(0.01)
    assert "b".upper() == "B"


def test_f():
    time.sleep(0.01)
    assert sorted([4, 1, 3, 2]) == [1, 2, 3, 4]
