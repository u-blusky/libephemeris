from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

import libephemeris.cli as cli_module
import libephemeris.dev_cli as dev_cli_module
import libephemeris.download as download_module
import libephemeris.iers_data as iers_data_module
import libephemeris.state as state_module


def _write_file(path: Path, size: int = 16) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * size)


def test_download_all_includes_iers_downloads(monkeypatch):
    calls: list[tuple] = []

    def fake_download_leb2_for_tier(
        tier_name: str,
        force: bool,
        show_progress: bool,
        quiet: bool,
        activate: bool,
    ) -> None:
        calls.append(("leb2", tier_name, force, show_progress, quiet, activate))

    def fake_download_for_tier(
        tier_name: str,
        force: bool,
        show_progress: bool,
        quiet: bool,
    ) -> None:
        calls.append(("tier", tier_name, force, show_progress, quiet))

    monkeypatch.setattr(
        download_module,
        "download_leb2_for_tier",
        fake_download_leb2_for_tier,
    )
    monkeypatch.setattr(download_module, "download_for_tier", fake_download_for_tier)
    monkeypatch.setattr(
        iers_data_module,
        "download_iers_finals",
        lambda force=False: calls.append(("iers_finals", force)),
    )
    monkeypatch.setattr(
        iers_data_module,
        "download_leap_seconds",
        lambda force=False: calls.append(("leap_seconds", force)),
    )
    monkeypatch.setattr(
        iers_data_module,
        "download_delta_t_data",
        lambda force=False: calls.append(("delta_t", force)),
    )

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        ["download", "all", "--force", "--no-progress", "--quiet"],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        ("leb2", "base", True, False, True, False),
        ("tier", "base", True, False, True),
        ("leb2", "medium", True, False, True, False),
        ("tier", "medium", True, False, True),
        ("leb2", "extended", True, False, True, False),
        ("tier", "extended", True, False, True),
        ("iers_finals", True),
        ("leap_seconds", True),
        ("delta_t", True),
    ]


def test_status_verbose_levels(monkeypatch, tmp_path):
    home_dir = tmp_path / "home"
    data_dir = home_dir / ".libephemeris"
    leb_dir = data_dir / "leb"
    spk_dir = data_dir / "spk"
    assist_dir = home_dir / ".libephemeris" / "assist"
    iers_dir = data_dir / "iers_cache"

    _write_file(data_dir / "de440s.bsp", 32)
    _write_file(data_dir / "planet_centers_base.bsp", 48)
    _write_file(leb_dir / "ephemeris_medium.leb", 64)
    _write_file(leb_dir / "base_core.leb", 12)
    _write_file(leb_dir / "base_asteroids.leb", 12)
    _write_file(leb_dir / "base_apogee.leb", 12)
    _write_file(leb_dir / "base_uranians.leb", 12)
    _write_file(spk_dir / "alpha.bsp", 24)
    _write_file(spk_dir / "beta.bsp", 24)
    _write_file(assist_dir / "linux_p1550p2650.440", 20)
    _write_file(assist_dir / "sb441-n16.bsp", 20)
    _write_file(iers_dir / "finals2000A.data", 10)
    _write_file(iers_dir / "leap_seconds.dat", 10)
    _write_file(iers_dir / "deltat.data", 10)

    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.setenv("LIBEPHEMERIS_LEB", str(leb_dir / "ephemeris_medium.leb"))
    monkeypatch.setattr(download_module, "get_data_dir", lambda: data_dir)
    monkeypatch.setattr(state_module, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(state_module, "get_precision_tier", lambda: "base")
    monkeypatch.setattr(state_module, "get_spk_cache_dir", lambda: str(spk_dir))
    monkeypatch.setattr(
        iers_data_module,
        "get_iers_cache_info",
        lambda: {
            "cache_dir": str(iers_dir),
            "finals_exists": True,
            "finals_age_days": 2.0,
            "leap_seconds_exists": True,
            "leap_seconds_age_days": 3.0,
            "delta_t_exists": True,
            "delta_t_age_days": 4.0,
        },
    )

    runner = CliRunner()

    base_result = runner.invoke(cli_module.cli, ["status"])
    assert base_result.exit_code == 0, base_result.output
    assert "Path:" not in base_result.output
    assert "leap_seconds.dat" in base_result.output

    verbose_result = runner.invoke(cli_module.cli, ["status", "-v"])
    assert verbose_result.exit_code == 0, verbose_result.output
    assert f"Path: {data_dir / 'de440s.bsp'}" in verbose_result.output
    assert f"Path: {leb_dir / 'base_core.leb'}" in verbose_result.output
    assert "base_core.leb" in verbose_result.output

    very_verbose_result = runner.invoke(cli_module.cli, ["status", "-vv"])
    assert very_verbose_result.exit_code == 0, very_verbose_result.output
    assert f"Path: {spk_dir / 'alpha.bsp'}" in very_verbose_result.output
    assert f"Path: {spk_dir / 'beta.bsp'}" in very_verbose_result.output

    json_result = runner.invoke(cli_module.cli, ["status", "--json"])
    assert json_result.exit_code == 0, json_result.output
    parsed = json.loads(json_result.output)
    assert parsed["data_directory"] == str(data_dir)
    assert parsed["calc_mode"] == "leb"
    assert parsed["precision_tier"] == "base"


def test_leph_status_delegates_to_shared_status_renderer(monkeypatch):
    calls: list[tuple[bool, int]] = []

    monkeypatch.setattr(
        download_module,
        "print_data_status",
        lambda *, as_json=False, verbose=0: calls.append((as_json, verbose)),
    )

    runner = CliRunner()
    result = runner.invoke(dev_cli_module.cli, ["status", "--json", "-vv"])

    assert result.exit_code == 0, result.output
    assert calls == [(True, 2)]
