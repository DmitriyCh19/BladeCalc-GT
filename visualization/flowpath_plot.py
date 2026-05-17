import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from dataclasses import replace

from visualization.flowpath_models import MachineFlowPathPlotData, FlowPathRow, FlowPathSection


def plot_machine_flowpath(data: MachineFlowPathPlotData) -> None:
    fig, (ax_flow, ax_T, ax_p) = plt.subplots(
        3, 1,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1, 1]}
    )
    _plot_channel(ax_flow, data.sections)
    _plot_rows(ax_flow, data.rows)
    _plot_stations(ax_T, ax_p, data)
    _set_flow_limits(ax_flow, data)

    ax_flow.set_title(data.name)
    ax_flow.set_ylabel("R, м")
    ax_flow.set_aspect("equal", adjustable="box")
    ax_flow.grid(True)

    ax_T.set_ylabel("T*, К")
    ax_T.grid(True)

    ax_p.set_ylabel("p*, Па")
    ax_p.set_xlabel("x, м")
    ax_p.grid(True)


    plt.tight_layout()
    plt.show()

def _row_polygon(row: FlowPathRow) -> list[tuple[float, float]]:
    return [
        (row.x0, row.hub_in / 2),
        (row.x1, row.hub_out / 2),
        (row.x1, row.tip_out / 2),
        (row.x0, row.tip_in / 2),
    ]

def _plot_rows(ax, rows: list[FlowPathRow]) -> None:
    for row in rows:
        poly = Polygon(_row_polygon(row), closed=True, fill=False, linewidth=1.8)
        ax.add_patch(poly)

        x_mid = 0.5 * (row.x0 + row.x1)
        y_mid = 0.25 * (row.hub_in + row.tip_in + row.hub_out + row.tip_out) / 2

        label = _row_label(row)
        ax.text(x_mid, y_mid, label, ha="center", va="center", fontsize=8)

    x_values = [x for row in rows for x in (row.x0, row.x1)]
    r_values = [d / 2 for row in rows for d in (row.hub_in, row.hub_out, row.tip_in, row.tip_out)]

    if x_values and r_values:
        ax.set_xlim(min(x_values), max(x_values))
        ax.set_ylim(min(r_values) * 0.95, max(r_values) * 1.05)


def _row_label(row: FlowPathRow) -> str:
    if row.row_type == "rotor":
        return f"РК {row.stage_index}"
    if row.row_type == "stator":
        return f"НА {row.stage_index}" if row.machine == "HPC" else f"СА {row.stage_index}"
    return f"{row.row_type} {row.stage_index}"

def _plot_stations(ax_T, ax_p, data: MachineFlowPathPlotData) -> None:
    x = [station.x for station in data.stations]
    T = [station.T_total for station in data.stations]
    p = [station.p_total for station in data.stations]

    ax_T.plot(x, T, marker="o", linewidth=1.8)
    ax_p.plot(x, p, marker="o", linewidth=1.8)

    for station in data.stations:
        ax_T.text(station.x, station.T_total, station.name, fontsize=7, rotation=45)
        ax_p.text(station.x, station.p_total, station.name, fontsize=7, rotation=45)

def _plot_channel(ax, sections: list[FlowPathSection]) -> None:
    if not sections:
        return

    x = [section.x for section in sections]
    hub = [section.hub / 2 for section in sections]
    tip = [section.tip / 2 for section in sections]

    ax.plot(x, hub, linewidth=2.2)
    ax.plot(x, tip, linewidth=2.2)
    ax.plot([x[0], x[0]], [hub[0], tip[0]], linewidth=2.2)
    ax.plot([x[-1], x[-1]], [hub[-1], tip[-1]], linewidth=2.2)

def _set_flow_limits(ax, data: MachineFlowPathPlotData) -> None:
    x_values = [section.x for section in data.sections]
    d_values = [d for section in data.sections for d in (section.hub, section.tip)]

    for row in data.rows:
        x_values.extend([row.x0, row.x1])
        d_values.extend([row.hub_in, row.hub_out, row.tip_in, row.tip_out])

    if x_values and d_values:
        r_values = [d / 2 for d in d_values]
        ax.set_xlim(min(x_values), max(x_values))
        ax.set_ylim(min(r_values) * 0.95, max(r_values) * 1.05)
