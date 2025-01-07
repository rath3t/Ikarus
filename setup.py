# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
# this file is modified from dumux https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/setup.py

try:
    from dune.packagemetadata import metaData
except ImportError:
    from packagemetadata import metaData
from skbuild import setup

# When building a new package, update the version numbers below and run:
# access  docker container with mounted repo to /tmp/Ikarus
# build _ikarus
# cd /tmp/Ikarus
# /dune/dune-common/build-cmake/run-in-dune-env pip install twine scikit-build
# git config --global --add safe.directory /tmp/Ikarus
# /dune/dune-common/build-cmake/run-in-dune-env python setup.py sdist
# /dune/dune-common/build-cmake/run-in-dune-env python -m twine upload dist/* --verbose

ikarusVersion = "0.4.1"
duneVersion = "2.11.0"

metadata = metaData(duneVersion)[1]
metadata["version"] = ikarusVersion

# refactor since ikarus pypi package already exists
metadata["name"] = "pyikarus"
setup(**metadata)
