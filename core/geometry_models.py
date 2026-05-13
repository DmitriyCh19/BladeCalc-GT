from dataclasses import dataclass

@dataclass
class SectionDiameters:
    hub: float
    mid: float
    tip: float


@dataclass
class MachineGeometry:
    inlet: SectionDiameters
    outlet: SectionDiameters

@dataclass
class VelocityTriangle:
    c: float
    c_a: float
    c_u: float

    w: float | None = None

    alpha_deg: float | None = None
    beta_deg: float | None = None

@dataclass
class RotorGeometry:
    inlet: SectionDiameters
    outlet: SectionDiameters

    blade_height_in: float
    blade_height_out: float

    blade_chord: float

    blade_count: int

    solidity: float