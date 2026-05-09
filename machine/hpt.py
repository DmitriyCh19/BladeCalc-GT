import math
from dataclasses import dataclass

from core.gas_func import q_lambda
from core.geometry import area_calc, diameter_from_area, calculate_section_diameters
from core.geometry_models import SectionDiameters, MachineGeometry
from configs.constants import *
from configs.modes import *

# =========================================================
# PARAMETERS
# =========================================================

@dataclass
class HPTParameters:
    mode: str

    # gas / thermodynamics
    T_gas: float
    p_gas: float

    T_out: float
    p_out: float

    T_cooling_air: float

    # stage
    z: int
    efficiency: float
    y: float
    theta: float

    # flow
    G_in: float
    G_cooling_rel: float
    G_hpc_in: float

    # blade
    density: float
    blade_force_coeff: float
    sigma_allow: float
    k_sigma: float

    # flow angles
    alpha_2_deg: float

    # reduced velocities
    lambda_2: float
    lambda_gas: float

    # work
    L_stage: float

class HighPressureTurbine:

    def __init__(self, params:HPTParameters):
        self.params = params

        self.geometry = None

        self.u_mid = None
        self.n_hp = None

        self.T_blade_cooling = None

        self.d_hub_out_rel  = None
        self.D_mean = None

    def calculate(self):
        c_ad = math.sqrt(2 * self.params.L_stage / self.params.efficiency)
        self.u_mid = self.params.y * c_ad * math.sqrt(1 / self.params.z)

        # температура лопатки ТВД
        T_blade = 0.95 * (self.params.T_out + self.u_mid ** 2 / (R_GAS * 2 * K_GAS / (K_GAS - 1))) # заменить комплекс на переменную теплоёмкость
        self.T_blade_cooling = T_blade / 0.95 - self.params.theta * (T_blade / 0.95 - self.params.T_cooling_air)

        # относительная высота лопатки ТВД
        blade_height_rel = (2 * (self.u_mid ** 2) * self.params.density * self.params.blade_force_coeff * self.params.k_sigma
                            ) / (self.params.sigma_allow* 1e6)
        
        q_lambda_2_hpt = q_lambda(lam=self.params.lambda_2, k=K_GAS) # заменить k
        alpha_2_hpt_rad = math.radians(self.params.alpha_2_deg)
        F_out_hpt = area_calc(G=(self.params.G_in + self.params.G_cooling_rel * self.params.G_hpc_in),
                               T=self.params.T_out, s=S_GAS, p=self.params.p_out, q=(q_lambda_2_hpt * math.sin(alpha_2_hpt_rad)))
        
        h_blabe_out_hpt = math.sqrt(F_out_hpt / (math.pi * blade_height_rel))

        hpt_out = calculate_section_diameters(
            D_ref=h_blabe_out_hpt * blade_height_rel,
            F=F_out_hpt,
            mode_name='mid',
            MODES_D=MODES_D
        )

        self.d_hub_out_rel = hpt_out.hub / hpt_out.tip

        q_lambda_g = q_lambda(lam=self.params.lambda_gas, k=K_GAS) # заменить k

        F_in_hpt = area_calc(G=self.params.G_in, T=self.params.T_gas, s=S_GAS, p=self.params.p_gas, q=q_lambda_g)

        # Базовые диаметры ТВД
        reference_diameters_hpt = {
            'hub': hpt_out.hub,
            'mid': hpt_out.mid,
            'tip': hpt_out.tip,
        }

        mode = MODES_D[self.params.mode]

        hpt_in = calculate_section_diameters(
            D_ref=reference_diameters_hpt[mode['ref']],
            F=F_in_hpt,
            mode_name=self.params.mode,
            MODES_D=MODES_D
        )

        self.n_hp = 60 * self.u_mid / (math.pi * hpt_out.mid)

        self.D_mean = (hpt_in.mid + hpt_out.mid) / 2

        self.geometry = MachineGeometry(
            inlet=hpt_in,
            outlet=hpt_out
        )