#!/usr/bin/env python3
"""
Version bumping script for fast-tcp package.
Updates version numbers in all relevant files.
"""

import argparse
import re
import sys
from pathlib import Path


def update_version_in_file(filepath, old_version, new_version):
    """Update version in a specific file."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()

        # Different patterns for different files
        patterns = [
            (r'__version__ = "[^"]*"', f'__version__ = "{new_version}"'),
            (r'version = "[^"]*"', f'version = "{new_version}"'),
            (r'version="[^"]*"', f'version="{new_version}"'),
        ]

        updated = False
        for pattern, replacement in patterns:
            if re.search(pattern, content):
                content = re.sub(pattern, replacement, content)
                updated = True

        if updated:
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(content)
            print(f"✓ Updated {filepath}")
            return True
        else:
            print(f"- No version found in {filepath}")
            return False

    except Exception as e:
        print(f"✗ Error updating {filepath}: {e}")
        return False


def get_current_version():
    """Get current version from __init__.py."""
    init_file = Path("fast_tcp/__init__.py")
    if not init_file.exists():
        print("Error: fast_tcp/__init__.py not found")
        return None

    with open(init_file, "r", encoding="utf-8") as f:
        content = f.read()

    match = re.search(r'__version__ = "([^"]*)"', content)
    if match:
        return match.group(1)

    print("Error: Could not find version in fast_tcp/__init__.py")
    return None


def parse_version(version_str):
    """Parse version string into components."""
    # Handle versions like 1.0.0, 1.0.0-rc.1, etc.
    match = re.match(r"^(\d+)\.(\d+)\.(\d+)(?:-(.+))?$", version_str)
    if not match:
        raise ValueError(f"Invalid version format: {version_str}")

    major, minor, patch, pre = match.groups()
    return int(major), int(minor), int(patch), pre


def bump_version(current_version, bump_type):
    """Bump version based on type."""
    major, minor, patch, pre = parse_version(current_version)

    if bump_type == "major":
        return f"{major + 1}.0.0"
    elif bump_type == "minor":
        return f"{major}.{minor + 1}.0"
    elif bump_type == "patch":
        return f"{major}.{minor}.{patch + 1}"
    elif bump_type == "rc":
        if pre and pre.startswith("rc."):
            rc_num = int(pre.split(".")[1]) + 1
            return f"{major}.{minor}.{patch}-rc.{rc_num}"
        else:
            return f"{major}.{minor}.{patch + 1}-rc.1"
    else:
        raise ValueError(f"Invalid bump type: {bump_type}")


def main():
    parser = argparse.ArgumentParser(
        description="Bump version for fast-tcp package",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/bump_version.py patch          # 1.0.0 -> 1.0.1
  python scripts/bump_version.py minor          # 1.0.0 -> 1.1.0
  python scripts/bump_version.py major          # 1.0.0 -> 2.0.0
  python scripts/bump_version.py rc             # 1.0.0 -> 1.0.1-rc.1
  python scripts/bump_version.py --version 1.2.3  # Set specific version
        """,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "bump_type",
        nargs="?",
        choices=["major", "minor", "patch", "rc"],
        help="Type of version bump",
    )
    group.add_argument("--version", type=str, help="Set specific version (e.g., 1.2.3)")

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be changed without making changes",
    )

    args = parser.parse_args()

    # Get current version
    current_version = get_current_version()
    if not current_version:
        sys.exit(1)

    print(f"Current version: {current_version}")

    # Calculate new version
    if args.version:
        new_version = args.version
        try:
            parse_version(new_version)  # Validate format
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        try:
            new_version = bump_version(current_version, args.bump_type)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)

    print(f"New version: {new_version}")

    if args.dry_run:
        print("\nDry run - no files will be changed")
        return

    # Update files
    files_to_update = [
        "fast_tcp/__init__.py",
        "setup.py",
        "pyproject.toml",
    ]

    print("\nUpdating files:")
    success_count = 0
    for filepath in files_to_update:
        if Path(filepath).exists():
            if update_version_in_file(filepath, current_version, new_version):
                success_count += 1
        else:
            print(f"- {filepath} not found")

    print(f"\nSuccessfully updated {success_count} files")

    if success_count > 0:
        print(f"\nNext steps:")
        print(f"1. git add .")
        print(f"2. git commit -m 'Bump version to {new_version}'")
        print(f"3. git tag v{new_version}")
        print(f"4. git push origin main")
        print(f"5. git push origin v{new_version}")


if __name__ == "__main__":
    main()
