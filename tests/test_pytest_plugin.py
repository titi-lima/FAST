from fast_tcp.integrations.pytest.plugin import pytest_collection_modifyitems


class _DummyConfig:
    def getoption(self, name: str):
        raise ValueError(f"no option named {name!r}")


class _DummySession:
    config = None


def test_collection_modifyitems_ignores_missing_fast_tcp_option() -> None:
    sentinel = object()
    items = [sentinel]
    pytest_collection_modifyitems(_DummySession(), _DummyConfig(), items)
    assert items == [sentinel]
