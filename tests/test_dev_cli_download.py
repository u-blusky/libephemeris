from __future__ import annotations

from click.testing import CliRunner

import libephemeris.dev_cli as dev_cli_module
import libephemeris.dev_cli.cmd_download as cmd_download_module


def test_leph_download_help_shows_dev_prerequisites() -> None:
    runner = CliRunner()
    result = runner.invoke(dev_cli_module.cli, ["download", "--help"])

    assert result.exit_code == 0, result.output
    assert "all" in result.output
    assert "planet-centers-sources" in result.output
    assert "iers" in result.output
    assert "leb-medium" not in result.output


def test_leph_download_all_runs_prerequisite_sequence(monkeypatch) -> None:
    calls: list[tuple] = []

    monkeypatch.setattr(
        cmd_download_module,
        "_download_spk_tier",
        lambda tier, force=False: calls.append(("spk", tier, force)),
    )
    monkeypatch.setattr(
        cmd_download_module,
        "_download_spk_extended",
        lambda force=False: calls.append(("spk_extended", force)),
    )
    monkeypatch.setattr(
        cmd_download_module,
        "_download_iers",
        lambda force=False: calls.append(("iers", force)),
    )
    monkeypatch.setattr(
        cmd_download_module,
        "_download_planet_center_sources",
        lambda force=False, tier=None: calls.append(("planet_centers", tier, force)),
    )
    monkeypatch.setattr(
        cmd_download_module,
        "_download_assist",
        lambda force=False: calls.append(("assist", force)),
    )

    runner = CliRunner()
    result = runner.invoke(dev_cli_module.cli, ["download", "all", "--force"])

    assert result.exit_code == 0, result.output
    assert calls == [
        ("spk", "base", True),
        ("spk", "medium", True),
        ("spk_extended", True),
        ("iers", True),
        ("planet_centers", None, True),
        ("assist", True),
    ]


def test_leph_download_planet_centers_sources_uses_download_only_mode(
    monkeypatch,
) -> None:
    calls: list[list[str]] = []

    monkeypatch.setattr(
        cmd_download_module,
        "_run_python",
        lambda args: calls.append(args) or 0,
    )

    runner = CliRunner()
    result = runner.invoke(
        dev_cli_module.cli,
        ["download", "planet-centers-sources", "--tier", "medium", "--force"],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        [
            "scripts/generate_planet_centers_spk.py",
            "--download-only",
            "--tier",
            "medium",
            "--force",
        ]
    ]
