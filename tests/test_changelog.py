"""
Tests for CHANGELOG.md validation.

Verifies that:
1. CHANGELOG.md exists and is properly formatted
2. Follows Keep a Changelog format
3. Contains required sections (Added, Changed, Fixed, etc.)
4. Version entries match pyproject.toml version
5. Links are properly formatted
"""

import os
import re
from datetime import datetime

import pytest


class TestChangelogExists:
    """Test that CHANGELOG.md exists and is readable."""

    def test_changelog_file_exists(self):
        """CHANGELOG.md should exist in project root."""
        changelog_path = os.path.join(os.path.dirname(__file__), "..", "CHANGELOG.md")
        assert os.path.exists(changelog_path), (
            f"CHANGELOG.md not found at {changelog_path}"
        )

    def test_changelog_is_not_empty(self):
        """CHANGELOG.md should have substantial content."""
        changelog_path = os.path.join(os.path.dirname(__file__), "..", "CHANGELOG.md")
        if os.path.exists(changelog_path):
            with open(changelog_path, "r", encoding="utf-8") as f:
                content = f.read()
            # Should have at least 1KB of content
            assert len(content) > 1000, (
                f"CHANGELOG.md seems too short: {len(content)} bytes"
            )


class TestChangelogFormat:
    """Test that CHANGELOG follows Keep a Changelog format."""

    @pytest.fixture
    def changelog_content(self):
        """Load CHANGELOG.md content."""
        changelog_path = os.path.join(os.path.dirname(__file__), "..", "CHANGELOG.md")
        with open(changelog_path, "r", encoding="utf-8") as f:
            return f.read()

    def test_has_changelog_header(self, changelog_content):
        """CHANGELOG should have a header."""
        assert "# Changelog" in changelog_content, (
            "CHANGELOG.md should have '# Changelog' header"
        )

    def test_has_keep_a_changelog_reference(self, changelog_content):
        """CHANGELOG should reference Keep a Changelog format."""
        assert "Keep a Changelog" in changelog_content, (
            "CHANGELOG.md should reference Keep a Changelog format"
        )

    def test_has_semver_reference(self, changelog_content):
        """CHANGELOG should reference Semantic Versioning."""
        assert "Semantic Versioning" in changelog_content, (
            "CHANGELOG.md should reference Semantic Versioning"
        )

    def test_has_unreleased_section(self, changelog_content):
        """CHANGELOG should have an Unreleased section."""
        assert "## [Unreleased]" in changelog_content, (
            "CHANGELOG.md should have an [Unreleased] section"
        )

    def test_has_version_entries(self, changelog_content):
        """CHANGELOG should have at least one version entry."""
        # Match version pattern like [0.1.8] - 2025-01-26
        version_pattern = r"## \[\d+\.\d+\.\d+\] - \d{4}-\d{2}-\d{2}"
        matches = re.findall(version_pattern, changelog_content)
        assert len(matches) >= 1, "CHANGELOG.md should have at least one version entry"

    def test_version_dates_are_valid(self, changelog_content):
        """Version dates should be valid ISO 8601 format (YYYY-MM-DD)."""
        # Find all dates in version headers
        date_pattern = r"## \[\d+\.\d+\.\d+\] - (\d{4}-\d{2}-\d{2})"
        matches = re.findall(date_pattern, changelog_content)

        for date_str in matches:
            try:
                datetime.strptime(date_str, "%Y-%m-%d")
            except ValueError:
                pytest.fail(f"Invalid date format in CHANGELOG: {date_str}")


class TestChangelogSections:
    """Test that CHANGELOG has proper section types."""

    @pytest.fixture
    def changelog_content(self):
        """Load CHANGELOG.md content."""
        changelog_path = os.path.join(os.path.dirname(__file__), "..", "CHANGELOG.md")
        with open(changelog_path, "r", encoding="utf-8") as f:
            return f.read()

    def test_has_added_section(self, changelog_content):
        """CHANGELOG should have an Added section."""
        assert "### Added" in changelog_content, (
            "CHANGELOG.md should have an 'Added' section"
        )

    def test_uses_valid_section_types(self, changelog_content):
        """All H3 sections should be valid Keep a Changelog types."""
        valid_types = [
            "Added",
            "Changed",
            "Deprecated",
            "Removed",
            "Fixed",
            "Security",
            "Documentation",
            "Tests",
        ]

        # Find all H3 headers (but not H4 sub-sections)
        h3_pattern = r"^### (\w+)"
        matches = re.findall(h3_pattern, changelog_content, re.MULTILINE)

        for section in matches:
            assert section in valid_types, (
                f"Invalid section type in CHANGELOG: '{section}'. "
                f"Valid types are: {valid_types}"
            )


class TestChangelogVersionSync:
    """Test that CHANGELOG version matches pyproject.toml."""

    def test_current_version_in_changelog(self):
        """Current version from pyproject.toml should be in CHANGELOG."""
        # Read pyproject.toml version
        pyproject_path = os.path.join(os.path.dirname(__file__), "..", "pyproject.toml")
        with open(pyproject_path, "r", encoding="utf-8") as f:
            pyproject_content = f.read()

        version_match = re.search(r'version = "(\d+\.\d+\.\d+)"', pyproject_content)
        assert version_match, "Could not find version in pyproject.toml"
        current_version = version_match.group(1)

        # Check CHANGELOG has this version
        changelog_path = os.path.join(os.path.dirname(__file__), "..", "CHANGELOG.md")
        with open(changelog_path, "r", encoding="utf-8") as f:
            changelog_content = f.read()

        assert f"[{current_version}]" in changelog_content, (
            f"Current version {current_version} not found in CHANGELOG.md"
        )


class TestChangelogLinks:
    """Test that CHANGELOG has proper comparison links."""

    @pytest.fixture
    def changelog_content(self):
        """Load CHANGELOG.md content."""
        changelog_path = os.path.join(os.path.dirname(__file__), "..", "CHANGELOG.md")
        with open(changelog_path, "r", encoding="utf-8") as f:
            return f.read()

    def test_has_unreleased_link(self, changelog_content):
        """CHANGELOG should have an Unreleased comparison link."""
        assert "[Unreleased]:" in changelog_content, (
            "CHANGELOG.md should have an [Unreleased] comparison link at the bottom"
        )

    def test_has_version_links(self, changelog_content):
        """CHANGELOG should have links for each version."""
        # Find all version headers
        version_pattern = r"## \[(\d+\.\d+\.\d+)\]"
        versions = re.findall(version_pattern, changelog_content)

        # Each version should have a link at the bottom
        for version in versions:
            link_pattern = f"[{version}]:"
            assert link_pattern in changelog_content, (
                f"Missing link for version {version} at bottom of CHANGELOG"
            )


class TestChangelogContent:
    """Test CHANGELOG content quality."""

    @pytest.fixture
    def changelog_content(self):
        """Load CHANGELOG.md content."""
        changelog_path = os.path.join(os.path.dirname(__file__), "..", "CHANGELOG.md")
        with open(changelog_path, "r", encoding="utf-8") as f:
            return f.read()

    def test_entries_are_bullet_points(self, changelog_content):
        """Changelog entries should use bullet points (-)."""
        # Check that there are bullet points in the content
        assert "\n- " in changelog_content, (
            "CHANGELOG entries should use bullet points (- item)"
        )

    def test_mentions_key_features(self, changelog_content):
        """CHANGELOG should mention key library features."""
        key_features = [
            "planetary",
            "house",
            "ayanamsha",
            "eclipse",
            "crossing",
        ]

        content_lower = changelog_content.lower()
        mentioned = [f for f in key_features if f in content_lower]

        assert len(mentioned) >= 3, (
            f"CHANGELOG should mention key features. Found only: {mentioned}"
        )

    def test_no_duplicate_versions(self, changelog_content):
        """Each version should only appear once as a header."""
        version_pattern = r"## \[(\d+\.\d+\.\d+)\]"
        versions = re.findall(version_pattern, changelog_content)

        seen = set()
        duplicates = []
        for version in versions:
            if version in seen:
                duplicates.append(version)
            seen.add(version)

        assert len(duplicates) == 0, f"Duplicate version entries found: {duplicates}"

    def test_versions_in_descending_order(self, changelog_content):
        """Version entries should be in descending order (newest first)."""
        version_pattern = r"## \[(\d+\.\d+\.\d+)\]"
        versions = re.findall(version_pattern, changelog_content)

        if len(versions) >= 2:
            # Parse versions as tuples for comparison
            def parse_version(v):
                return tuple(int(x) for x in v.split("."))

            for i in range(len(versions) - 1):
                current = parse_version(versions[i])
                next_v = parse_version(versions[i + 1])
                assert current >= next_v, (
                    f"Versions not in descending order: "
                    f"{versions[i]} should come after {versions[i + 1]}"
                )
