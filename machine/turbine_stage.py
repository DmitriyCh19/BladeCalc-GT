import math
from dataclasses import dataclass

from configs.constants import *
from configs.modes import *
from core.geometry_models import MachineGeometry

@dataclass
class TurbineStageParameters:
    mode: str

    G_in: float
    G_out: float
    h_2_rotor: float
    u_mid: float
    turbine: MachineGeometry

    p_in: float
    pi: float
    T_in: float
    T_out:float

    reaction: float
    y_1: float
    S_stator: float
    S_rotor: float

    r_stator: float
    r_rotor: float
    r_stator_colling: float | None = None
    r_rotor_colling: float | None = None

    phi_cooling: float | None = None
    psi_cooling: float | None = None
