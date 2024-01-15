#!/usr/bin/env python

# When specifying the path to be added to `sys.path` in a testing context,
# it's often more robust to use an absolute path rather than a relative path.
# This helps ensure that the correct directory is added, regardless of the 
# current working directory or the location from which the tests are run.
# `os.path.abspath` is used to obtain the absolute path. This approach is 
# especially useful in scenarios where the test runner might change the 
# current working directory before executing the tests. Using an absolute 
# path reduces the dependency on the current working directory and provides
# a more reliable way to construct the path to the project's parent
# directory, which is then added to `sys.path`.

# Add the parent directory to the sys.path
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from moldyn.base import *
