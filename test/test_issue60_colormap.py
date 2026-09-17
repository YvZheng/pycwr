"""Issue #60: bundled radar colormaps work without importing Py-ART."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest


_RADAR_SETUP = """
import numpy as np
from pycwr.core.NRadar import PRD
radar = PRD(
    fields={name: np.ones((4, 3)) for name in ("dBZ", "ZDR", "CC", "KDP", "PhiDP", "W")},
    scan_type="ppi", time=np.arange(4), range=np.array([1000., 2000., 3000.]),
    azimuth=np.array([0., 90., 180., 270.]), elevation=np.ones(4),
    latitude=30., longitude=120., altitude=100.,
    sweep_start_ray_index=np.array([0]), sweep_end_ray_index=np.array([3]),
    fixed_angle=np.array([1.]), bins_per_sweep=np.array([3]),
    nyquist_velocity=np.array([15.]), frequency=5.6,
    unambiguous_range=np.array([150000.]), nrays=4, nsweeps=1, sitename="TEST",
)
"""


class Issue60ColormapTests(unittest.TestCase):
    def _run_fresh(self, script):
        import pycwr

        package_file = Path(pycwr.__file__).resolve()
        env = os.environ.copy()
        env["PYTHONPATH"] = str(package_file.parent.parent)
        prelude = 'import matplotlib\nmatplotlib.use("Agg")\nimport matplotlib.pyplot as plt\n'
        prelude += (
            "import pycwr\nfrom pathlib import Path\n"
            f"assert Path(pycwr.__file__).resolve() == Path({str(package_file)!r})\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [sys.executable, "-c", prelude + textwrap.dedent(script)],
                cwd=tmpdir, env=env, capture_output=True, text=True, timeout=60,
            )
        self.assertEqual(result.returncode, 0, msg=result.stdout + result.stderr)

    def test_easy_ppi_dualpol_fields_without_pyart(self):
        self._run_fresh(_RADAR_SETUP + """
from pycwr.draw import plot_ppi
for field in ("ZDR", "CC", "KDP", "PhiDP", "W"):
    result = plot_ppi(radar, field=field, show=False)
    result.fig.canvas.draw()
    plt.close(result.fig)
assert "pyart" not in __import__("sys").modules
""")

    def test_direct_graph_ppi_zdr_without_pyart(self):
        self._run_fresh(_RADAR_SETUP + """
from pycwr.draw.RadarPlot import Graph
fig, ax = plt.subplots()
Graph(radar).plot_ppi(ax, 0, "ZDR")
fig.canvas.draw()
plt.close(fig)
assert "pyart" not in __import__("sys").modules
""")

    def test_all_default_colormaps_resolve_in_fresh_registry(self):
        self._run_fresh("""
import numpy as np
from pycwr.configure.default_config import CINRAD_COLORMAP
from pycwr.draw._plot_core import _resolve_matplotlib_cmap
for name in set(CINRAD_COLORMAP.values()):
    cmap = _resolve_matplotlib_cmap(name)
    if name.startswith("pyart_"):
        expected = plt.get_cmap("copy_" + name)
        np.testing.assert_array_equal(cmap(np.linspace(0, 1, 17)), expected(np.linspace(0, 1, 17)))
for name in ("pyart_RefDiff_r", "pyart_HomeyerRainbow", "pyart_HomeyerRainbow_r"):
    assert _resolve_matplotlib_cmap(name) is not None
assert "pyart" not in __import__("sys").modules
""")

    def test_existing_external_registration_is_preserved(self):
        self._run_fresh("""
import numpy as np
from matplotlib.colors import ListedColormap
external = ListedColormap(["red", "green"], name="pyart_RefDiff")
if hasattr(matplotlib, "colormaps"):
    matplotlib.colormaps.register(external)
else:
    matplotlib.cm.register_cmap(name="pyart_RefDiff", cmap=external)
from pycwr.draw import colormap
from pycwr.draw._plot_core import _resolve_matplotlib_cmap
np.testing.assert_array_equal(_resolve_matplotlib_cmap("pyart_RefDiff")([0, 1]), external([0, 1]))
np.testing.assert_array_equal(plt.get_cmap("pyart_RefDiff")([0, 1]), external([0, 1]))
""")

    def test_unknown_colormap_still_raises(self):
        self._run_fresh("""
from pycwr.draw._plot_core import _resolve_matplotlib_cmap
try:
    _resolve_matplotlib_cmap("pyart_not_a_real_colormap")
except ValueError:
    pass
else:
    raise AssertionError("unknown colormap must raise ValueError")
""")

    @unittest.skipUnless(importlib.util.find_spec("pyart") is not None, "external Py-ART is not installed")
    def test_pyart_import_order_does_not_break_zdr(self):
        for order in ("import pyart\nfrom pycwr.draw import colormap\n",
                      "from pycwr.draw import colormap\nimport pyart\n"):
            with self.subTest(order=order):
                self._run_fresh(order + _RADAR_SETUP + """
from pycwr.draw import plot_ppi
result = plot_ppi(radar, field="ZDR", show=False)
result.fig.canvas.draw()
plt.close(result.fig)
""")


if __name__ == "__main__":
    unittest.main()
