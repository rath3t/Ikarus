#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: CC0-1.0
import os
from datetime import datetime


def read_old_version():
    script_dir = os.path.dirname(os.path.realpath("__file__"))
    rel_path = "dune.module"
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
    version_number_array = version_number.split(".")

    if len(version_number_array) == 2:
        major, minor = version_number_array
        patch = 0
    elif len(version_number_array) == 3:
        major, minor, patch = version_number_array
    else:
        raise Exception("Bad versoin number passed!")
    return f"{major}.{minor}.{int(patch) + 1}"


def inplace_change(filename: str, old_string: str, new_string: str):
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            print('"{old_string}" not found in {filename}.'.format(**locals()))
            return

    with open(filename, "w") as f:
        print(
            'Changing "{old_string}" to "{new_string}" in {filename}'.format(**locals())
        )
        s = s.replace(old_string, new_string)
        f.write(s)


def changeLine(filename: str, old_string: str, new_string: str):
    with open(filename) as f:
        s = f.read()
        for line in s.split("\n"):
            print(line)
            if line.startswith(old_string):
                with open(filename, "w") as fw:
                    print(
                        'Changing the line "{line}"\n to "{new_string}" in {filename}'.format(
                            **locals()
                        )
                    )
                    s = s.replace(line, new_string)
                    fw.write(s)


def update_all_versions(version_override=None):
    """Update all version numbers in local files"""
    old_version_number = read_old_version()
    if version_override == "dev":
        new_version_number = bump_patch_number(old_version_number)
    elif version_override is None:
        new_version_number = old_version_number
    else:
        new_version_number = version_override

    if version_override == "dev":
        new_version_number += ".dev" + datetime.now().strftime("%Y%m%d%H%M%S")

    print(f"Bump version from {old_version_number} to {new_version_number}")
    if version_override != "dev":
        inplace_change(
            "dune.module",
            f"Version: {old_version_number}",
            f"Version: {new_version_number}",
        )
        inplace_change(
            "CMakeLists.txt",
            f"VERSION {old_version_number}",
            f"VERSION {new_version_number}",
        )

    changeLine(
        "setup.py",
        f"ikarusVersion =",
        f'ikarusVersion = "{new_version_number}"',
    )


import sys

if __name__ == "__main__":
    try:
        var = sys.argv[1]
        var = var.removeprefix("v")
        update_all_versions(var)
    except IndexError:
        update_all_versions()
