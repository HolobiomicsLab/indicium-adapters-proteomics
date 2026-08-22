"""Tests that the adapter is registered and discoverable."""

from __future__ import annotations

from importlib.metadata import entry_points

from indicium_adapters import CURRENT_API_VERSION, discover_adapters


def test_entry_point_registered():
    names = {
        ep.name
        for ep in entry_points(group="indicium.adapters")
        if ep.module.startswith("indicium_adapters_proteomics")
    }
    assert names == {"proteomics"}


def test_discover_adapters_finds_proteomics():
    found = discover_adapters()
    assert "proteomics" in found, f"discover_adapters() missing 'proteomics'; found={sorted(found)}"


def test_discovered_adapter_has_expected_api_version():
    assert discover_adapters()["proteomics"].api_version == CURRENT_API_VERSION
