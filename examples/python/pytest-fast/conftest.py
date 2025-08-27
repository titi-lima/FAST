import pytest


def pytest_configure(config: "pytest.Config") -> None:
    """Ensure the installed FAST TCP plugin is available without duplicating options.

    If the entry-point plugin is already loaded (name: 'fast_tcp'), do nothing.
    Otherwise, register the in-package plugin module explicitly.
    """
    pm = config.pluginmanager
    if pm.hasplugin("fast_tcp"):
        return
    try:
        from fast_tcp.integrations.pytest import plugin as fast_tcp_plugin
    except Exception:
        return
    pm.register(fast_tcp_plugin, name="fast_tcp")
