import math
from dataclasses import dataclass

from core.gas_func import q_lambda
from core.geometry import area_calc, calculate_section_diameters
from core.geometry_models import MachineGeometry
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

        self.blade_height_rel = None
        self.h_blade_out = None
        self.d_hub_out_rel  = None
        self.D_mean = None

    def calculate(self):
        c_ad = math.sqrt(2 * self.params.L_stage / self.params.efficiency)
        self.u_mid = self.params.y * c_ad * math.sqrt(1 / self.params.z)

        # температура лопатки ТВД
        T_blade = 0.95 * (self.params.T_out + self.u_mid ** 2 / (R_GAS * 2 * K_GAS / (K_GAS - 1))) # заменить комплекс на переменную теплоёмкость
        self.T_blade_cooling = T_blade / 0.95 - self.params.theta * (T_blade / 0.95 - self.params.T_cooling_air)
        
        q_lambda_2 = q_lambda(lam=self.params.lambda_2, k=K_GAS) # заменить k
        alpha_2_rad = math.radians(self.params.alpha_2_deg)
        F_out = area_calc(G=(self.params.G_in + self.params.G_cooling_rel * self.params.G_hpc_in),
                               T=self.params.T_out, s=S_GAS, p=self.params.p_out, q=(q_lambda_2 * math.sin(alpha_2_rad)))
        
        if self.params.z == 1:
            self.blade_height_out_rel = (2 * (self.u_mid ** 2) * self.params.density * self.params.blade_force_coeff * self.params.k_sigma
                            ) / (self.params.sigma_allow* 1e6)
            self.h_blade_out = math.sqrt(F_out / (math.pi * self.blade_height_out_rel))
        elif self.params.z == 2:
            sigma_r_1 = 0.85 * self.params.sigma_allow / self.params.k_sigma
            k_sigma_1 = self.params.sigma_allow / sigma_r_1
            if k_sigma_1 >= 2 and k_sigma_1 <= 1.8:
                ValueError(f"recomend change params: 1.8 <= {k_sigma_1} <= 2.0")
            self.blade_height_out_rel = (2 * (self.u_mid ** 2) * self.params.density * self.params.blade_force_coeff * self.params.k_sigma
                            ) / (self.params.sigma_allow* 1e6)
            self.h_blade_out = math.sqrt(F_out / (math.pi * self.blade_height_out_rel))
            self.blade_height_1_rel = (2 * (self.u_mid ** 2) * self.params.density * self.params.blade_force_coeff * k_sigma_1
                            ) / (self.params.sigma_allow* 1e6)
            self.h_blade_1 = self.h_blade_out * 0.85
            self.F_1_out = self.h_blade_1 * math.pi * self.blade_height_1_rel
        else:  
            raise ValueError(f"{self.params.z} > 2 WIP")

        hpt_out = calculate_section_diameters(
            D_ref=self.h_blade_out * self.blade_height_out_rel,
            F=F_out,
            mode_name='mid',
            MODES_D=MODES_D
        )

        self.d_hub_out_rel = hpt_out.hub / hpt_out.tip

        q_lambda_g = q_lambda(lam=self.params.lambda_gas, k=K_GAS) # заменить k

        F_in = area_calc(G=self.params.G_in, T=self.params.T_gas, s=S_GAS, p=self.params.p_gas, q=q_lambda_g)

        # Базовые диаметры ТВД
        reference_diameters_hpt = {
            'hub': hpt_out.hub,
            'mid': hpt_out.mid,
            'tip': hpt_out.tip,
        }

        mode = MODES_D[self.params.mode]

        hpt_in = calculate_section_diameters(
            D_ref=reference_diameters_hpt[mode['ref']],
            F=F_in,
            mode_name=self.params.mode,
            MODES_D=MODES_D
        )

        self.n_hp = 60 * self.u_mid / (math.pi * hpt_out.mid)

        self.D_mean = (hpt_in.mid + hpt_out.mid) / 2

        self.geometry = MachineGeometry(
            inlet=hpt_in,
            outlet=hpt_out
        )
