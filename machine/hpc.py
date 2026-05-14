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


class HighPressureCompressor:

    def __init__(self, params:HPCParameters):
        self.params = params

        self.geometry = None

        self.u_tip_in = None
        self.c_out = None
        self.h_blade_out = None

        self.D_mean = None

        self.stage_machines = []
        self.stage_results = []
        self.stages_result = None

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
    def calculate_stages(
            self,
            *,
            alpha1_deg: float,
            L_stage_rel: list[float],
            eff_stage: list[float],
            eff_rel: list[float],
            reaction_stage: list[float],
            c1a_stage: list[float],
            h_rot_rel: list[float],
            p_in: float,
            D_const: float | None = None,
            D_tip: float | None = None,
            G: float | None = None,
            d_hub_rel: float | None = None,
            lambda_in: float | None = None,
            u_tip: float | None = None,
            T_in: float | None = None,
            D_1mid: float | None = None,
            D_1hub: float | None = None,
            K_g: float | None = None,
            eff_rotor: float = 0.93,
            mode: str | None = None,
            z: int | None = None,
    ) -> HPCStagesResult:
        """
        Последовательный расчёт всех ступеней компрессора.

        На вход первой ступени подаются T_in и p_in.
        Для каждой следующей ступени:
            T_in = T_out предыдущей ступени,
            p_in = p_out предыдущей ступени.

        После расчёта рабочих колёс выполняется расчёт направляющих аппаратов
        между соседними ступенями.
        """

        if z is None:
            z = len(L_stage_rel)

        self._check_stage_input_lengths(
            z=z,
            L_stage_rel=L_stage_rel,
            eff_stage=eff_stage,
            eff_rel=eff_rel,
            reaction_stage=reaction_stage,
            c1a_stage=c1a_stage,
            h_rot_rel=h_rot_rel,
        )

        mode = self.params.mode if mode is None else mode
        G = self.params.G_in if G is None else G
        d_hub_rel = self.params.d_hub_in_rel if d_hub_rel is None else d_hub_rel
        lambda_in = self.params.lambda_in if lambda_in is None else lambda_in
        K_g = self.params.K_g if K_g is None else K_g
        T_current = self.params.T_in if T_in is None else T_in
        p_current = p_in

        if D_const is None:
            if self.geometry is None:
                raise ValueError(
                    "D_const не задан, а self.geometry ещё не рассчитана. "
                    "Сначала нужно вызвать hpc.calculate(n) или явно передать D_const."
                )
            D_const = getattr(self.geometry.inlet, mode)

        if D_tip is None:
            D_tip = D_const

        if u_tip is None:
            if self.u_tip_in is None:
                raise ValueError(
                    "u_tip не задан, а self.u_tip_in ещё не рассчитан. "
                    "Сначала нужно вызвать hpc.calculate(n) или явно передать u_tip."
                )
            u_tip = self.u_tip_in

        if D_1mid is None:
            D_1mid = self.geometry.inlet.mid if self.geometry is not None else 0.0

        if D_1hub is None:
            D_1hub = self.geometry.inlet.hub if self.geometry is not None else 0.0

        self.stage_machines = []
        self.stage_results = []

        for i in range(z):
            params_stage = CompressorStageParameters(
                mode=mode,
                D_const=D_const,

                alpha1_deg=alpha1_deg,
                L_rel=L_stage_rel[i],
                D_tip=D_tip,
                efficiency=eff_stage[i],
                eff_tip_rel=eff_rel[i],
                reaction=reaction_stage[i],
                c_1a=c1a_stage[i],
                c_3a=c1a_stage[i + 1],
                h_rot_rel=h_rot_rel[i],
                G=G,
                d_hub_rel=d_hub_rel,
                lambda_in=lambda_in,
                u_tip=u_tip,
                T_in=T_current,
                p_in=p_current,
                D_1mid=D_1mid,
                D_1hub=D_1hub,
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

        self.stages_result = HPCStagesResult(
            stages=self.stage_results,
            T_out=self.stage_results[-1].thermodynamics.T_out,
            p_out=self.stage_results[-1].thermodynamics.p_out,
            pi_total=pi_total,
            L_total=L_total,
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
                raise ValueError(
                    f"Массив {name} содержит {len(values)} элементов, "
                    f"а требуется минимум {z}."
                )

        if len(c1a_stage) < z + 1:
            raise ValueError(
                f"Массив c1a_stage содержит {len(c1a_stage)} элементов, "
                f"а требуется минимум {z + 1}, потому что для каждой ступени "
                f"используются c1a_stage[i] и c1a_stage[i + 1]."
            )


