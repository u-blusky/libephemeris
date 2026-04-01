from __future__ import annotations

from click.testing import CliRunner

import libephemeris.dev_cli as dev_cli_module
import libephemeris.dev_cli.cmd_download as cmd_download_module


def test_leph_download_help_lists_all_command() -> None:
    runner = CliRunner()
    result = runner.invoke(dev_cli_module.cli, ["download", "--help"])

    assert result.exit_code == 0, result.output
    assert "all" in result.output
    assert "Full developer dataset" in result.output


def test_leph_download_all_runs_full_developer_sequence(monkeypatch) -> None:
    calls: list[list[str]] = []

    monkeypatch.setattr(
        cmd_download_module,
        "_run_python",
        lambda args: calls.append(args) or 0,
    )

    runner = CliRunner()
    result = runner.invoke(
        dev_cli_module.cli,
        ["download", "all", "--force", "--no-progress", "--quiet"],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        [
            "-m",
            "libephemeris.cli",
            "download",
            "all",
            "--force",
            "--no-progress",
            "--quiet",
        ],
        [
            "-m",
            "libephemeris.cli",
            "download",
            "leb-base",
            "--force",
            "--no-progress",
            "--quiet",
        ],
        [
            "-m",
            "libephemeris.cli",
            "download",
            "leb-medium",
            "--force",
            "--no-progress",
            "--quiet",
        ],
        [
            "-m",
            "libephemeris.cli",
            "download",
            "leb-extended",
            "--force",
            "--no-progress",
            "--quiet",
        ],
        [
            "-m",
            "libephemeris.cli",
            "download",
            "assist",
            "--force",
            "--no-progress",
            "--quiet",
        ],
    ]
