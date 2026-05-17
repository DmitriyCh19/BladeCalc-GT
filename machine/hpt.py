import math
from dataclasses import dataclass

from core.gas_func import q_lambda
from core.geometry import area_calc, calculate_section_diameters
from core.geometry_models import MachineGeometry
from configs.constants import *
from configs.modes import *
from machine.turbine_stage import TurbineStage
from machine.turbine_stage_models import (
    TurbineStageParameters, TurbineStageFlowParams, TurbineStageThermoParams,
    TurbineStageKinematicParams, TurbineStageGeometryParams, TurbineStageBladeRowsParams,
    BladeRowInput, TurbineStageResult
)

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

@dataclass
class HPTStagesResult:
    stages: list[TurbineStageResult]
    T_out: float
    p_out: float
    pi_total: float
    L_total: float
    length_total: float
    checks: dict[str, float]

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

        self.G_out = None
        self.F_out = None
        self.stage_machines = []
        self.stage_results = []
        self.stages_result = None

    def calculate(self):
        c_ad = math.sqrt(2 * self.params.L_stage / self.params.efficiency)
        self.u_mid = self.params.y * c_ad * math.sqrt(1 / self.params.z)

        # температура лопатки ТВД
        T_blade = 0.95 * (self.params.T_out + self.u_mid ** 2 / (R_GAS * 2 * K_GAS / (K_GAS - 1))) # заменить комплекс на переменную теплоёмкость
        self.T_blade_cooling = T_blade / 0.95 - self.params.theta * (T_blade / 0.95 - self.params.T_cooling_air)
        
        q_lambda_2 = q_lambda(lam=self.params.lambda_2, k=K_GAS) # заменить k
        alpha_2_rad = math.radians(self.params.alpha_2_deg)

        self.G_out = self.params.G_in + self.params.G_cooling_rel * self.params.G_hpc_in
        F_out = area_calc(G=self.G_out, T=self.params.T_out, s=S_GAS,
        p=self.params.p_out, q=q_lambda_2 * math.sin(alpha_2_rad))
        self.F_out = F_out
        
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
    def calculate_stages(
        self, *,
        L_stage: list[float],
        eff_stage: list[float],
        reaction_stage: list[float],
        stator_axial_chord: list[float],
        rotor_axial_chord: list[float],
        stator_trailing_edge_radius_rel: list[float],
        rotor_trailing_edge_radius_rel: list[float],
        cooling_rel_stage: list[float],
        y_stage: list[float] | None = None,
        pi_stage: list[float] | None = None,
        phi_stage: list[float] | None = None,
        psi_stage: list[float] | None = None,
        stator_cooling_pitch_factor: list[float] | None = None,
        rotor_cooling_pitch_factor: list[float] | None = None,
    ) -> HPTStagesResult:
        z = self.params.z

        if self.geometry is None or self.u_mid is None or self.G_out is None:
            raise ValueError("Перед calculate_stages() нужно вызвать hpt.calculate().")
        if z > 2:
            raise ValueError(f"{z} > 2 WIP")

        h_2_rotor = self._get_stage_rotor_out_heights()
        h_2_rotor_rel = self._get_stage_rotor_out_height_rels()

        self._check_stage_input_lengths(
            z=z, L_stage=L_stage, eff_stage=eff_stage, reaction_stage=reaction_stage,
            stator_axial_chord=stator_axial_chord, rotor_axial_chord=rotor_axial_chord,
            stator_trailing_edge_radius_rel=stator_trailing_edge_radius_rel,
            rotor_trailing_edge_radius_rel=rotor_trailing_edge_radius_rel,
            cooling_rel_stage=cooling_rel_stage,
        )

        self._check_cooling_sum(cooling_rel_stage=cooling_rel_stage)

        y_stage = self._stage_list(y_stage, z, self.params.y)
        phi_stage = self._stage_list(phi_stage, z, 1.0)
        psi_stage = self._stage_list(psi_stage, z, 1.0)
        stator_cooling_pitch_factor = self._stage_list(stator_cooling_pitch_factor, z, 1.0)
        rotor_cooling_pitch_factor = self._stage_list(rotor_cooling_pitch_factor, z, 1.0)

        if pi_stage is None:
            pi_stage = [(self.params.p_gas / self.params.p_out) ** (1 / z) for _ in range(z)]

        T_current = self.params.T_gas
        p_current = self.params.p_gas
        G_current = self.params.G_in
        inlet_current = self.geometry.inlet

        self.stage_machines = []
        self.stage_results = []

        for i in range(z):
            G_out_current = G_current + cooling_rel_stage[i] * self.params.G_hpc_in
            T_out_target = self._stage_T_out(T_in=T_current, L_stage=L_stage[i])

            params_stage = TurbineStageParameters(
                flow=TurbineStageFlowParams(G_in=G_current, G_out=G_out_current),
                thermo=TurbineStageThermoParams(
                    p_in=p_current, pi=pi_stage[i], T_in=T_current,
                    T_out=T_out_target, L_stage=L_stage[i],
                    efficiency=eff_stage[i]
                ),
                kinematics=TurbineStageKinematicParams(
                    u_mid=self.u_mid, reaction=reaction_stage[i],
                    phi_cooling=phi_stage[i], psi_cooling=psi_stage[i]
                ),
                geometry=TurbineStageGeometryParams(
                    mode=self.params.mode,inlet=inlet_current,
                    h_2_rotor=h_2_rotor[i],h_2_rotor_rel=h_2_rotor_rel[i],
                ),
                blade_rows=TurbineStageBladeRowsParams(
                    stator=BladeRowInput(
                        axial_chord=stator_axial_chord[i],
                        trailing_edge_radius_rel=stator_trailing_edge_radius_rel[i],
                        cooling_pitch_factor=stator_cooling_pitch_factor[i],
                    ),
                    rotor=BladeRowInput(
                        axial_chord=rotor_axial_chord[i],
                        trailing_edge_radius_rel=rotor_trailing_edge_radius_rel[i],
                        cooling_pitch_factor=rotor_cooling_pitch_factor[i],
                    ),
                ),
            )

            stage_machine = TurbineStage(params_stage)
            stage_result = stage_machine.calculate()

            self.stage_machines.append(stage_machine)
            self.stage_results.append(stage_result)

            T_current = stage_result.thermodynamics.T_2
            p_current = stage_result.thermodynamics.p_out
            G_current = G_out_current
            inlet_current = stage_result.geometry.rotor_outlet

        length_total = sum(stage.length.total for stage in self.stage_results)
        checks = self.check_stage_outlet_match()

        self.stages_result = HPTStagesResult(
            stages=self.stage_results,
            T_out=self.stage_results[-1].thermodynamics.T_2,
            p_out=self.stage_results[-1].thermodynamics.p_out,
            pi_total=self.params.p_gas / self.stage_results[-1].thermodynamics.p_out,
            L_total=sum(stage.thermodynamics.L_stage for stage in self.stage_results),
            length_total=length_total, checks=checks,
        )

        return self.stages_result
    
    def _get_stage_rotor_out_height_rels(self) -> list[float]:
        if self.params.z == 1:
            return [self.blade_height_out_rel]
        if self.params.z == 2:
            return [self.blade_height_1_rel, self.blade_height_out_rel]
        raise ValueError(f"{self.params.z} > 2 WIP")
    
    def _get_stage_rotor_out_heights(self) -> list[float]:
        if self.params.z == 1:
            return [self.h_blade_out]
        if self.params.z == 2:
            return [self.h_blade_1, self.h_blade_out]
        raise ValueError(f"{self.params.z} > 2 WIP")

    @staticmethod
    def _stage_list(values: list[float] | None, z: int, default: float) -> list[float]:
        if values is None:
            return [default for _ in range(z)]
        if len(values) < z:
            raise ValueError(f"Массив содержит {len(values)} элементов, а требуется минимум {z}.")
        return values

    @staticmethod
    def _stage_T_out(T_in: float, L_stage: float) -> float:
        cp_gas = K_GAS / (K_GAS - 1) * R_GAS
        return T_in - L_stage / cp_gas

    def _check_cooling_sum(self, cooling_rel_stage: list[float]) -> None:
        cooling_sum = sum(cooling_rel_stage)
        if abs(cooling_sum - self.params.G_cooling_rel) > 1e-6:
            raise ValueError(
                f"Сумма cooling_rel_stage должна быть равна G_cooling_rel: "
                f"{cooling_sum:.6f} != {self.params.G_cooling_rel:.6f}"
            )

    @staticmethod
    def _check_stage_input_lengths(
        *, z: int, L_stage: list[float], reaction_stage: list[float],
        stator_axial_chord: list[float], rotor_axial_chord: list[float],
        stator_trailing_edge_radius_rel: list[float], eff_stage: list[float],
        rotor_trailing_edge_radius_rel: list[float],
        cooling_rel_stage: list[float],
    ) -> None:
        arrays = {
            "L_stage": L_stage,
            "eff_stage": eff_stage,
            "reaction_stage": reaction_stage,
            "stator_axial_chord": stator_axial_chord,
            "rotor_axial_chord": rotor_axial_chord,
            "stator_trailing_edge_radius_rel": stator_trailing_edge_radius_rel,
            "rotor_trailing_edge_radius_rel": rotor_trailing_edge_radius_rel,
            "cooling_rel_stage": cooling_rel_stage,
        }

        for name, values in arrays.items():
            if len(values) < z:
                raise ValueError(f"Массив {name} содержит {len(values)} элементов, а требуется минимум {z}.")
            

    def check_stage_outlet_match(self) -> dict[str, float]:
        if not self.stage_results:
            raise ValueError("Нет рассчитанных ступеней ТВД.")

        def rel_error_percent(calc: float, ref: float) -> float:
            if abs(ref) < 1e-12:
                return 0.0 if abs(calc) < 1e-12 else float("inf")
            return (calc - ref) / ref * 100

        last = self.stage_results[-1]
        G_stage_out = self.params.G_in + self.params.G_cooling_rel * self.params.G_hpc_in

        return {
            "dG_out_%": rel_error_percent(G_stage_out, self.G_out),
            "dT_out_%": rel_error_percent(last.thermodynamics.T_2, self.params.T_out),
            "dp_out_%": rel_error_percent(last.thermodynamics.p_out, self.params.p_out),
            "dalpha_out_%": rel_error_percent(last.velocity.rotor_outlet.alpha_deg, self.params.alpha_2_deg),
            "dlambda_out_%": rel_error_percent(last.velocity.lambda_c2, self.params.lambda_2),
            "dF_out_%": rel_error_percent(last.geometry.F_2, self.F_out),
            "dD_hub_out_%": rel_error_percent(last.geometry.rotor_outlet.hub, self.geometry.outlet.hub),
            "dD_mid_out_%": rel_error_percent(last.geometry.rotor_outlet.mid, self.geometry.outlet.mid),
            "dD_tip_out_%": rel_error_percent(last.geometry.rotor_outlet.tip, self.geometry.outlet.tip),
        }