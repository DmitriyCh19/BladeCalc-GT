import math
from dataclasses import dataclass

from core.geometry_models import SectionDiameters, MachineGeometry
from core.gas_func import q_lambda, polotrop, velocity_critical
from core.geometry import area_calc, diameter_from_area, reconstruct_section_diameters, calculate_section_diameters
from configs.constants import *
from configs.modes import *

@dataclass
class HPCParameters:
    mode: str

    # gas
    T_in: float
    
    T_out: float
    p_out: float

    pi: float

    d_hub_in_rel: float
    K_g: float

    efficiency: float
    lambda_in : float
    lambda_out: float

    # flow
    G_in : float
    G_cooling_rel: float


class HighPressureCompressor:

    def __init__(self, params:HPCParameters):
        self.params = params

        self.geometry = None

        self.u_tip_in = None
        self.c_out = None
        self.h_blade_out = None

        self.D_mean = None

    def calculate(self, n:float):
        a_crit_out = velocity_critical(k=K_AIR, R=R_AIR, T=self.params.T_out)
        self.c_out = self.params.lambda_out * a_crit_out

        q_lambda_out_hpc = q_lambda(lam=self.params.lambda_out, k=K_AIR)
        q_lambda_in_hpc = q_lambda(lam=self.params.lambda_in, k=K_AIR)
        F_out = area_calc(G=(self.params.G_in - self.params.G_cooling_rel * self.params.G_in), T=self.params.T_out, s=S_AIR, p=self.params.p_out, q=(q_lambda_out_hpc * self.params.K_g))

        politropa = polotrop(pi=self.params.pi, T_in=self.params.T_in, T_out=self.params.T_out)

        F_rel = self.params.pi ** ((politropa + 1) / (2 * politropa)) * q_lambda_out_hpc / q_lambda_in_hpc
        F_in = F_rel * F_out

        hpc_in = reconstruct_section_diameters(F=F_in, hub_to_tip_ratio=self.params.d_hub_in_rel)

        reference_diameters_hpc = {
            'hub': hpc_in.hub,
            'mid': hpc_in.mid,
            'tip': hpc_in.tip,
        }

        mode = MODES_D[self.params.mode]

        hpc_out = calculate_section_diameters(
            D_ref=reference_diameters_hpc[mode['ref']],
            F=F_out,
            mode_name=self.params.mode,
            MODES_D=MODES_D
        )

        self.h_blade_out = (hpc_out.tip - hpc_out.hub) / 2

        self.u_tip_in = math.pi * hpc_in.tip * n / 60

        self.D_mean = (hpc_in.mid + hpc_out.mid) / 2

        self.geometry = MachineGeometry(
            inlet=hpc_in,
            outlet=hpc_out
        )


