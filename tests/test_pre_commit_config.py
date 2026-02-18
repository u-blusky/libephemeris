"""Tests for pre-commit configuration validation.

These tests verify that the pre-commit configuration file is valid
and correctly configured with all required hooks.

Note: These tests are skipped because .pre-commit-config.yaml is not
maintained in this project.
"""

import pytest

pytestmark = pytest.mark.skip(reason="pre-commit config not maintained in this project")

from pathlib import Path

import yaml


@pytest.fixture
def pre_commit_config_path():
    """Get the path to the pre-commit config file."""
    return Path(__file__).parent.parent / ".pre-commit-config.yaml"


@pytest.fixture
def pre_commit_config(pre_commit_config_path):
    """Load and parse the pre-commit configuration."""
    with open(pre_commit_config_path) as f:
        return yaml.safe_load(f)


class TestPreCommitConfigExists:
    """Test that the pre-commit configuration file exists."""

    def test_pre_commit_config_file_exists(self, pre_commit_config_path):
        """Verify .pre-commit-config.yaml exists in the project root."""
        assert pre_commit_config_path.exists(), (
            f"Pre-commit config file not found at {pre_commit_config_path}"
        )


class TestPreCommitConfigStructure:
    """Test the structure of the pre-commit configuration."""

    def test_config_has_repos_key(self, pre_commit_config):
        """Verify the config contains the required 'repos' key."""
        assert "repos" in pre_commit_config, "Config must contain 'repos' key"

    def test_repos_is_list(self, pre_commit_config):
        """Verify 'repos' is a list."""
        assert isinstance(pre_commit_config["repos"], list), "'repos' must be a list"

    def test_repos_not_empty(self, pre_commit_config):
        """Verify there are repositories configured."""
        assert len(pre_commit_config["repos"]) > 0, (
            "At least one repo must be configured"
        )


class TestRequiredHooks:
    """Test that all required hooks are configured."""

    def get_all_hook_ids(self, pre_commit_config):
        """Extract all hook IDs from the configuration."""
        hook_ids = []
        for repo in pre_commit_config["repos"]:
            if "hooks" in repo:
                for hook in repo["hooks"]:
                    if "id" in hook:
                        hook_ids.append(hook["id"])
        return hook_ids

    def test_ruff_hook_configured(self, pre_commit_config):
        """Verify ruff linting hook is configured."""
        hook_ids = self.get_all_hook_ids(pre_commit_config)
        assert "ruff" in hook_ids, "Ruff linting hook must be configured"

    def test_ruff_format_hook_configured(self, pre_commit_config):
        """Verify ruff-format hook is configured."""
        hook_ids = self.get_all_hook_ids(pre_commit_config)
        assert "ruff-format" in hook_ids, "Ruff format hook must be configured"

    def test_black_hook_configured(self, pre_commit_config):
        """Verify black formatting hook is configured."""
        hook_ids = self.get_all_hook_ids(pre_commit_config)
        assert "black" in hook_ids, "Black formatting hook must be configured"

    def test_mypy_hook_configured(self, pre_commit_config):
        """Verify mypy type checking hook is configured."""
        hook_ids = self.get_all_hook_ids(pre_commit_config)
        assert "mypy" in hook_ids, "Mypy type checking hook must be configured"

    def test_pytest_hook_configured(self, pre_commit_config):
        """Verify pytest hook is configured."""
        hook_ids = self.get_all_hook_ids(pre_commit_config)
        assert "pytest" in hook_ids, "Pytest hook must be configured"


class TestRuffHookConfiguration:
    """Test ruff hook specific configuration."""

    def get_ruff_repo(self, pre_commit_config):
        """Get the ruff repository configuration."""
        for repo in pre_commit_config["repos"]:
            if "repo" in repo and "ruff" in repo["repo"]:
                return repo
        return None

    def test_ruff_repo_exists(self, pre_commit_config):
        """Verify ruff repository is configured."""
        repo = self.get_ruff_repo(pre_commit_config)
        assert repo is not None, "Ruff repository must be configured"

    def test_ruff_has_version(self, pre_commit_config):
        """Verify ruff repository has a version pinned."""
        repo = self.get_ruff_repo(pre_commit_config)
        assert repo is not None, "Ruff repo must exist"
        assert "rev" in repo, "Ruff repo must have 'rev' (version) specified"

    def test_ruff_has_fix_arg(self, pre_commit_config):
        """Verify ruff hook is configured with --fix argument."""
        repo = self.get_ruff_repo(pre_commit_config)
        assert repo is not None, "Ruff repo must exist"
        ruff_hook = None
        for hook in repo["hooks"]:
            if hook["id"] == "ruff":
                ruff_hook = hook
                break
        assert ruff_hook is not None, "Ruff hook not found"
        assert "args" in ruff_hook, "Ruff hook should have args"
        assert "--fix" in ruff_hook["args"], "Ruff hook should include --fix"


class TestBlackHookConfiguration:
    """Test black hook specific configuration."""

    def get_black_repo(self, pre_commit_config):
        """Get the black repository configuration."""
        for repo in pre_commit_config["repos"]:
            if "repo" in repo and "black" in repo["repo"]:
                return repo
        return None

    def test_black_repo_exists(self, pre_commit_config):
        """Verify black repository is configured."""
        repo = self.get_black_repo(pre_commit_config)
        assert repo is not None, "Black repository must be configured"

    def test_black_has_version(self, pre_commit_config):
        """Verify black repository has a version pinned."""
        repo = self.get_black_repo(pre_commit_config)
        assert repo is not None, "Black repo must exist"
        assert "rev" in repo, "Black repo must have 'rev' (version) specified"


class TestMypyHookConfiguration:
    """Test mypy hook specific configuration."""

    def get_mypy_repo(self, pre_commit_config):
        """Get the mypy repository configuration."""
        for repo in pre_commit_config["repos"]:
            if "repo" in repo and "mypy" in repo["repo"]:
                return repo
        return None

    def test_mypy_repo_exists(self, pre_commit_config):
        """Verify mypy repository is configured."""
        repo = self.get_mypy_repo(pre_commit_config)
        assert repo is not None, "Mypy repository must be configured"

    def test_mypy_has_version(self, pre_commit_config):
        """Verify mypy repository has a version pinned."""
        repo = self.get_mypy_repo(pre_commit_config)
        assert repo is not None, "Mypy repo must exist"
        assert "rev" in repo, "Mypy repo must have 'rev' (version) specified"


class TestPytestHookConfiguration:
    """Test pytest hook specific configuration."""

    def get_pytest_hook(self, pre_commit_config):
        """Get the pytest hook configuration."""
        for repo in pre_commit_config["repos"]:
            if "hooks" in repo:
                for hook in repo["hooks"]:
                    if hook.get("id") == "pytest":
                        return hook
        return None

    def test_pytest_hook_exists(self, pre_commit_config):
        """Verify pytest hook is configured."""
        hook = self.get_pytest_hook(pre_commit_config)
        assert hook is not None, "Pytest hook must be configured"

    def test_pytest_is_local_hook(self, pre_commit_config):
        """Verify pytest is configured as a local hook."""
        for repo in pre_commit_config["repos"]:
            if repo.get("repo") == "local":
                for hook in repo.get("hooks", []):
                    if hook.get("id") == "pytest":
                        return  # Found it in local repo
        pytest.fail("Pytest hook should be in a local repository")

    def test_pytest_excludes_slow_tests(self, pre_commit_config):
        """Verify pytest hook excludes slow tests for fast pre-commit."""
        hook = self.get_pytest_hook(pre_commit_config)
        assert hook is not None, "Pytest hook must exist"
        entry = hook.get("entry", "")
        assert "not slow" in entry or "-m" in entry, (
            "Pytest hook should exclude slow tests for fast pre-commit"
        )
