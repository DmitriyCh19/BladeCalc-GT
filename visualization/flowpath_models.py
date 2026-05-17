from dataclasses import dataclass

@dataclass
class FlowPathRow:
    machine: str
    stage_index: int
    row_type: str      # "rotor", "stator", "nozzle"
    x0: float
    x1: float
    hub_in: float
    hub_out: float
    tip_in: float
    tip_out: float

@dataclass
class FlowPathStation:
    name: str
    x: float
    T_total: float | None
    p_total: float | None

@dataclass
class MachineFlowPathPlotData:
    name: str
    sections: list[FlowPathSection]
    rows: list[FlowPathRow]
    stations: list[FlowPathStation]
    x_max: float

@dataclass
class FlowPathSection:
    x: float
    hub: float
    tip: float