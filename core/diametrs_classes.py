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