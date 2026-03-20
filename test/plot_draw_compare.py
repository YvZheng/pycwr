import argparse
from pathlib import Path
import sys
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

warnings.filterwarnings(
    "ignore",
    message="The input coordinates to pcolormesh are interpreted as cell centers.*",
)

from pycwr.core.transforms import antenna_vectors_to_cartesian_cwr, cartesian_to_geographic_aeqd
from pycwr.draw._plot_core import (
    CartesianReferenceOptions,
    ColorbarOptions,
    MapOptions,
    PlotStyle,
    ccrs,
    plot_cartesian,
    plot_map,
    plot_vertical_section,
)


def legacy_fix_ticks(ticks):
    if (ticks % 1).sum() == 0:
        return ["%2.f" % i for i in ticks]
    return ["%.2f" % i for i in ticks]


def legacy_set_grids(ax, max_range_km, delta_r_km=50.0):
    theta = np.linspace(0, 2 * np.pi, 200)
    rings = np.arange(delta_r_km, max_range_km + 1, delta_r_km)
    for radius in rings:
        x0 = radius * np.cos(theta)
        y0 = radius * np.sin(theta)
        ax.plot(x0, y0, linestyle="-", linewidth=0.6, color="#5B5B5B")
    for rad in np.arange(0, np.pi, np.pi / 6):
        ax.plot(
            [-1 * rings[-1] * np.sin(rad), rings[-1] * np.sin(rad)],
            [-1 * rings[-1] * np.cos(rad), rings[-1] * np.cos(rad)],
            linestyle="-",
            linewidth=0.6,
            color="#5B5B5B",
        )


def legacy_plot_ppi(ax, x_m, y_m, radar_data, title, vmin=-5.0, vmax=75.0, cmap="CN_ref", cmap_bins=16):
    levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
    norm = BoundaryNorm(levels, ncolors=plt.get_cmap(cmap).N, clip=True)
    mesh = ax.pcolormesh(x_m / 1000.0, y_m / 1000.0, radar_data, cmap=cmap, norm=norm, shading="auto")
    range_cycle = float(np.max(x_m) / 1000.0)
    legacy_set_grids(ax, range_cycle)
    ax.set_aspect("equal")
    ax.set_xlim([-range_cycle, range_cycle])
    ax.set_ylim([-range_cycle, range_cycle])
    ax.set_xlabel("Distance From Radar In East (Uints:km)")
    ax.set_ylabel("Distance From Radar In North (Uints:km)")
    ax.set_title(title)
    cbar = ax.figure.colorbar(mesh, ax=ax, ticks=levels)
    cbar.set_ticklabels(legacy_fix_ticks(levels))
    cbar.set_label("reflectivity (dBZ)")
    return mesh


def legacy_plot_section(ax, mesh_xy_m, mesh_z_m, field_data, title, height_km=(0, 18), vmin=-5.0, vmax=75.0, cmap="CN_ref", cmap_bins=16):
    levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
    norm = BoundaryNorm(levels, ncolors=plt.get_cmap(cmap).N, clip=True)
    mesh = ax.pcolormesh(
        mesh_xy_m / 1000.0,
        mesh_z_m / 1000.0,
        field_data[:, :-1],
        cmap=cmap,
        norm=norm,
        shading="auto",
    )
    ax.set_xlabel("Distance From Section Start (Uints:km)")
    ax.set_ylabel("Height (Uints:km)")
    ax.set_title(title)
    ax.set_ylim(list(height_km))
    ax.grid(True, zorder=15)
    cbar = ax.figure.colorbar(mesh, ax=ax, orientation="horizontal", ticks=levels, pad=0.16)
    cbar.set_ticklabels(legacy_fix_ticks(levels))
    cbar.set_label("reflectivity")
    return mesh


def legacy_plot_map(ax, lon, lat, radar_data, title, vmin=-5.0, vmax=75.0, cmap="CN_ref", cmap_bins=16):
    levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
    norm = BoundaryNorm(levels, ncolors=plt.get_cmap(cmap).N, clip=True)
    mesh = ax.pcolormesh(
        lon,
        lat,
        radar_data,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        norm=norm,
        zorder=4,
        shading="auto",
    )
    ax.set_extent([float(np.min(lon)), float(np.max(lon)), float(np.min(lat)), float(np.max(lat))], ccrs.PlateCarree())
    ax.set_title(title)
    cbar = ax.figure.colorbar(mesh, ax=ax, ticks=levels)
    cbar.set_ticklabels(legacy_fix_ticks(levels))
    cbar.set_label("reflectivity (dBZ)")
    return mesh


def build_demo_fields():
    ranges = np.linspace(1000.0, 150000.0, 160)
    azimuth = np.linspace(0.0, 359.0, 360)
    elevation = np.full_like(azimuth, 1.2)
    x_m, y_m, _ = antenna_vectors_to_cartesian_cwr(ranges, azimuth, elevation, h=0.0)
    radius_km = np.hypot(x_m / 1000.0, y_m / 1000.0)
    angle = np.deg2rad(np.broadcast_to(azimuth[:, None], x_m.shape))
    radar_data = (
        45.0 * np.exp(-((radius_km - 55.0) ** 2) / (2.0 * 18.0 ** 2))
        + 12.0 * np.cos(3.0 * angle)
        + 4.0 * np.sin(radius_km / 8.0)
    )
    radar_data = np.clip(radar_data, -5.0, 75.0)

    section_distance_m = np.linspace(0.0, 180000.0, 220)
    section_height_m = np.linspace(0.0, 18000.0, 120)
    mesh_xy_m = np.broadcast_to(section_distance_m[None, :], (section_height_m.size, section_distance_m.size))
    mesh_z_m = np.broadcast_to(section_height_m[:, None], (section_height_m.size, section_distance_m.size))
    section_field = (
        50.0 * np.exp(-((mesh_xy_m / 1000.0 - 85.0) ** 2) / (2.0 * 24.0 ** 2))
        * np.exp(-((mesh_z_m / 1000.0 - 6.5) ** 2) / (2.0 * 2.3 ** 2))
    )
    section_field = np.pad(section_field, ((0, 0), (0, 1)), mode="edge")

    station_lonlat = (118.78, 32.04)
    lon, lat = cartesian_to_geographic_aeqd(x_m, y_m, station_lonlat[0], station_lonlat[1])
    return x_m, y_m, radar_data, mesh_xy_m, mesh_z_m, section_field, lon, lat


def build_comparison_figure(output_path, include_map=False):
    x_m, y_m, radar_data, mesh_xy_m, mesh_z_m, section_field, lon, lat = build_demo_fields()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    legacy_plot_ppi(axes[0, 0], x_m, y_m, radar_data, "Legacy PPI")
    modern_ppi_style = PlotStyle(
        cmap="CN_ref",
        value_range=(-5.0, 75.0),
        bins=16,
        colorbar=ColorbarOptions(label="reflectivity (dBZ)"),
    )
    plot_cartesian(
        axes[0, 1],
        x_m,
        y_m,
        radar_data,
        modern_ppi_style,
        fig=fig,
        title="Refactored PPI",
        reference_options=CartesianReferenceOptions(ring_spacing_km=25.0, ring_color="#444444", spoke_color="#444444"),
    )

    legacy_plot_section(axes[1, 0], mesh_xy_m, mesh_z_m, section_field, "Legacy Vertical Section")
    modern_section_style = PlotStyle(
        cmap="CN_ref",
        value_range=(-5.0, 75.0),
        bins=16,
        colorbar=ColorbarOptions(orientation="horizontal", label="reflectivity (dBZ)"),
    )
    plot_vertical_section(
        axes[1, 1],
        [mesh_xy_m],
        [mesh_z_m],
        [section_field],
        modern_section_style,
        fig=fig,
        title="Refactored Vertical Section",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    if include_map and ccrs is not None:
        map_output = output_path.with_name(output_path.stem + "_map" + output_path.suffix)
        map_fig = plt.figure(figsize=(14, 5), constrained_layout=True)
        legacy_ax = map_fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        modern_ax = map_fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        legacy_plot_map(legacy_ax, lon, lat, radar_data, "Legacy Map PPI")
        modern_map_style = PlotStyle(
            cmap="CN_ref",
            value_range=(-5.0, 75.0),
            bins=16,
            colorbar=ColorbarOptions(label="reflectivity (dBZ)"),
        )
        plot_map(
            modern_ax,
            lon,
            lat,
            radar_data,
            modern_map_style,
            fig=map_fig,
            map_options=MapOptions(data_crs=ccrs.PlateCarree(), tick_step_degrees=0.5),
        )
        modern_ax.set_title("Refactored Map PPI")
        map_fig.savefig(map_output, dpi=150)
        plt.close(map_fig)
        return output_path, map_output
    return output_path, None


def main():
    parser = argparse.ArgumentParser(description="Generate before/after draw-module comparison figures.")
    parser.add_argument(
        "--output",
        default="test/artifacts/draw_compare.png",
        help="Output image path for the Cartesian and section comparison figure.",
    )
    parser.add_argument(
        "--include-map",
        action="store_true",
        help="Also generate a map comparison figure. This may require cartopy background data.",
    )
    args = parser.parse_args()
    output_path, map_output = build_comparison_figure(Path(args.output), include_map=args.include_map)
    print(output_path)
    if map_output is not None:
        print(map_output)


if __name__ == "__main__":
    main()
