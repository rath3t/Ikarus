#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: CC0-1.0
import os


def read_old_version():
    script_dir = os.path.dirname(os.path.realpath("__file__"))
    rel_path = "../../../dune.module"
    abs_file_path = os.path.join(script_dir, rel_path)
    print(abs_file_path)
    with open(abs_file_path) as f:
        s = f.read()
        for line in s.split("\n"):
            print(line)
            if line.startswith("Version: "):
                _, old_version = line.split(" ")
                return old_version


def bump_patch_number(version_number: str) -> str:
    """Return a copy of `version_number` with the patch number incremented."""
    major, minor, patch = version_number.split(".")
    return f"{major}.{minor}.{int(patch) + 1}"


def inplace_change(filename: str, old_string: str, new_string: str):
    # Safely read the input filename using 'with'
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            print('"{old_string}" not found in {filename}.'.format(**locals()))
            return

    # Safely write the changed content, if found in the file
    with open(filename, "w") as f:
        print(
            'Changing "{old_string}" to "{new_string}" in {filename}'.format(**locals())
        )
        s = s.replace(old_string, new_string)
        f.write(s)


def update_all_versions(version_override=None):
    """Update all version numbers in local files"""
    old_version_number = read_old_version()
    if version_override is None:
        new_version_number = bump_patch_number(old_version_number)
    else:
        new_version_number = version_override
    print(f"Bump version from {old_version_number} to {new_version_number}")
    inplace_change(
        "../../../dune.module",
        f"Version: {old_version_number}",
        f"Version: {new_version_number}",
    )
    inplace_change(
        "../../../CMakeLists.txt",
        f"VERSION {old_version_number}",
        f"VERSION {new_version_number}",
    )
    inplace_change(
        "../../../setup.py",
        f'ikarusVersion = "{old_version_number}"',
        f'ikarusVersion = "{new_version_number}"',
    )


import sys

if __name__ == "__main__":
    try:
        var = sys.argv[1]
        update_all_versions(var)
    except IndexError:
        update_all_versions()