# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import os
import re
import datetime


def update_license_header(file_path):
    try:
        with open(file_path, "r", newline="") as file:
            content = file.read()

        # Update the license header using regex
        updated_content, nsubs = re.subn(
            r"SPDX-FileCopyrightText: (\d{4}-\d{4}) The Ikarus Developers",
            lambda match: f"SPDX-FileCopyrightText: 2021-{datetime.datetime.now().year} The Ikarus Developers",
            content,
        )
        if nsubs > 0:
            with open(file_path, "w", newline="") as file:
                file.write(updated_content)
            print(f"Updated license header in: {file_path}")
        else:
            print(f"Didn't update license header in: {file_path}")

    except Exception as e:
        print(f"Error updating license header in {file_path}: {e}")


def update_dep5(file_path):
    try:
        with open(file_path, "r", newline="") as file:
            content = file.read()

        updated_content, nsubs = re.subn(
            r"Copyright: (\d{4}-\d{4}) The Ikarus Developers",
            lambda match: f"Copyright: 2021-{datetime.datetime.now().year} The Ikarus Developers",
            content,
        )
        if nsubs > 0:
            with open(file_path, "w", newline="") as file:
                file.write(updated_content)

    except Exception as e:
        print(f"Error updating license information in {file_path}: {e}")


def process_directory(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if os.path.basename(file_path) == "dep5":
                update_dep5(file_path)
            else:
                update_license_header(file_path)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath("__file__"))

    print(script_dir)
    process_directory(script_dir)
