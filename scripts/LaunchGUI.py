# -*- coding: utf-8 -*-
"""Launch the lightweight local web viewer."""

import os
import sys
import warnings

warnings.filterwarnings("ignore")

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from pycwr.GraphicalInterface.web_app import launch
except ImportError as exc:
    raise ImportError(
        "Flask is required to launch the web viewer. "
        "Install it with `pip install \"pycwr[full]\"` or `pip install Flask`."
    ) from exc


if __name__ == "__main__":
    launch()
