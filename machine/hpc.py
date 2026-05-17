import math
from dataclasses import dataclass

from core.geometry_models import MachineGeometry
from core.gas_func import q_lambda, polotrop, velocity_critical
from core.geometry import area_calc, reconstruct_section_diameters, calculate_section_diameters
from configs.constants import *
from configs.modes import *
from machine.compressor_stage import CompressorStageParameters, CompressorStage, CompressorStageResult

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

@dataclass
class HPCStagesResult:
    stages: list[CompressorStageResult]

    T_out: float
    p_out: float

    pi_total: float
    L_total: float
    checks: dict[str, float] | None = None


class HighPressureCompressor:

    def __init__(self, params:HPCParameters):
        self.params = params

        self.geometry = None

        self.u_tip_in = None
        self.c_out = None
        self.h_blade_out = None

        self.D_mean = None
        self.F_out = None

        self.stage_machines = []
        self.stage_results = []
        self.stages_result = None

    def calculate(self, n:float):
        a_crit_out = velocity_critical(k=K_AIR, R=R_AIR, T=self.params.T_out)
        self.c_out = self.params.lambda_out * a_crit_out

        q_lambda_out_hpc = q_lambda(lam=self.params.lambda_out, k=K_AIR)
        q_lambda_in_hpc = q_lambda(lam=self.params.lambda_in, k=K_AIR)

        G_out = self.params.G_in * (1 - self.params.G_cooling_rel)
        F_out = area_calc(G=G_out, T=self.params.T_out, s=S_AIR, p=self.params.p_out, q=q_lambda_out_hpc * self.params.K_g)
        self.F_out = F_out

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

    def calculate_stages(
        self, *,
        L_stage_rel: list[float],
        eff_stage: list[float],
        eff_rel: list[float],
        reaction_stage: list[float],
        c1a_stage: list[float],
        h_rot_rel: list[float],
        eff_rotor: float = 0.93,
        z: int | None = None,
    ) -> HPCStagesResult:
        if z is None:
            z = len(L_stage_rel)

        if self.geometry is None or self.u_tip_in is None:
            raise ValueError("Перед calculate_stages() нужно вызвать hpc.calculate(n).")

        self._check_stage_input_lengths(
            z=z, L_stage_rel=L_stage_rel, eff_stage=eff_stage,
            eff_rel=eff_rel, reaction_stage=reaction_stage,
            c1a_stage=c1a_stage, h_rot_rel=h_rot_rel,
        )

        mode = self.params.mode
        D_const = getattr(self.geometry.inlet, mode)
        G = self.params.G_in
        u_tip = self.u_tip_in
        T_current = self.params.T_in
        p_current = self.params.p_out / self.params.pi
        K_g = self.params.K_g

        self.stage_machines = []
        self.stage_results = []

        for i in range(z):
            params_stage = CompressorStageParameters(
                mode=mode,
                D_const=D_const,
                L_rel=L_stage_rel[i],
                efficiency=eff_stage[i],
                eff_tip_rel=eff_rel[i],
                reaction=reaction_stage[i],
                c_1a=c1a_stage[i],
                c_3a=c1a_stage[i + 1],
                h_rot_rel=h_rot_rel[i],
                G=G,
                u_tip=u_tip,
                T_in=T_current,
                p_in=p_current,
                K_g=K_g,
                eff_rotor=eff_rotor,
            )

            stage_machine = CompressorStage(params_stage)
            stage_machine.calculate()

            self.stage_machines.append(stage_machine)
            self.stage_results.append(stage_machine.stage)

            T_current = stage_machine.stage.thermodynamics.T_out
            p_current = stage_machine.stage.thermodynamics.p_out

        self.calculate_all_stators()

        pi_total = 1.0
        L_total = 0.0

        for stage in self.stage_results:
            pi_total *= stage.thermodynamics.pi
            L_total += stage.thermodynamics.L_stage

        checks = self.check_stage_outlet_match()

        self.stages_result = HPCStagesResult(
            stages=self.stage_results,
            T_out=self.stage_results[-1].thermodynamics.T_out,
            p_out=self.stage_results[-1].thermodynamics.p_out,
            pi_total=pi_total,
            L_total=L_total,
            checks=checks,
        )

        return self.stages_result
    
    def calculate_all_stators(self) -> None:
        """
        Расчёт НА после каждой ступени компрессора.

        Для всех ступеней, кроме последней, выходной угол НА принимается равным
        входному углу следующей ступени.

        Для последней ступени выходной угол НА принимается равным 90 градусов.
        """

        z = len(self.stage_machines)

        for i in range(z):
            current_stage = self.stage_machines[i]

            if i < z - 1:
                next_stage = self.stage_machines[i + 1].stage

                current_stage.calculate_stator(
                    alpha_out_deg=next_stage.inlet_triangle.alpha_deg
                )
            else:
                current_stage.calculate_stator(
                    alpha_out_deg=90.0
                )

    @staticmethod
    def _check_stage_input_lengths(
        *,
        z: int,
        L_stage_rel: list[float],
        eff_stage: list[float],
        eff_rel: list[float],
        reaction_stage: list[float],
        c1a_stage: list[float],
        h_rot_rel: list[float],
    ) -> None:
        arrays = {
            "L_stage_rel": L_stage_rel,
            "eff_stage": eff_stage,
            "eff_rel": eff_rel,
            "reaction_stage": reaction_stage,
            "h_rot_rel": h_rot_rel,
        }

        for name, values in arrays.items():
            if len(values) < z:
                raise ValueError(f"Массив {name} содержит {len(values)} элементов, а требуется минимум {z}.")

        if len(c1a_stage) < z + 1:
            raise ValueError(
                f"Массив c1a_stage содержит {len(c1a_stage)} элементов, "
                f"а требуется минимум {z + 1}, потому что используются c1a_stage[i] и c1a_stage[i + 1]."
            )

    def check_stage_outlet_match(self) -> dict[str, float]:
        if not self.stage_results:
            raise ValueError("Нет рассчитанных ступеней КВД.")
        if self.geometry is None or self.F_out is None:
            raise ValueError("Перед проверкой нужно вызвать hpc.calculate(n).")

        def rel_error_percent(calc: float, ref: float) -> float:
            if abs(ref) < 1e-12:
                return 0.0 if abs(calc) < 1e-12 else float("inf")
            return (calc - ref) / ref * 100

        last = self.stage_results[-1]

        return {
            "dT_out_%": rel_error_percent(last.thermodynamics.T_out, self.params.T_out),
            "dp_out_%": rel_error_percent(last.thermodynamics.p_out, self.params.p_out),
            "dF_out_%": rel_error_percent(last.F_outlet, self.F_out),
            "dh_blade_out_%": rel_error_percent(last.rotor.blade_height_out, self.h_blade_out),
            "dD_hub_out_%": rel_error_percent(last.rotor.outlet.hub, self.geometry.outlet.hub),
            "dD_mid_out_%": rel_error_percent(last.rotor.outlet.mid, self.geometry.outlet.mid),
            "dD_tip_out_%": rel_error_percent(last.rotor.outlet.tip, self.geometry.outlet.tip),
        }