import math
from dataclasses import dataclass

from configs.constants import *
from configs.modes import *
from machine.turbine_stage_models import (
    TurbineStageParameters, TurbineStageFlowParams, TurbineStageThermoParams,
    TurbineStageKinematicParams, TurbineStageGeometryParams, TurbineStageBladeRowsParams,
    BladeRowInput, TurbineStageWork, TurbineStatorOutlet, TurbineRotorOutlet,
    BladeRowGridParams, TurbineStageBladeRows, TurbineStageThermodynamics,
    TurbineStageGeometry, TurbineStageLoading, TurbineStageVelocity, TurbineStageResult
)
from core.gas_func import velocity_critical, pi_lambda, q_lambda
from core.geometry import area_calc, calculate_section_diameters
from core.geometry_models import VelocityTriangle



class TurbineStage:
    def __init__(self, params:TurbineStageParameters):
        self.params = params

    def calculate(self) -> TurbineStageResult:
        flow = self.params.flow
        thermo = self.params.thermo
        kin = self.params.kinematics
        geom = self.params.geometry
        rows = self.params.blade_rows

        work = self.calculate_stage_work(thermo=thermo, kin=kin)
        stator = self.calculate_stator_outlet(flow=flow, thermo=thermo, kin=kin, geom=geom, work=work)
        rotor = self.calculate_rotor_outlet(flow=flow, thermo=thermo, kin=kin, geom=geom, work=work, stator=stator)
        blade_rows = self.calculate_blade_rows(geom=geom, rows=rows, kin=kin, stator=stator, rotor=rotor)

        result = self.assemble_result(thermo=thermo, work=work, stator=stator, rotor=rotor, blade_rows=blade_rows)
        self.result = result
        return result

    def calculate_stage_work(self, thermo: TurbineStageThermoParams, kin: TurbineStageKinematicParams) -> TurbineStageWork:
        if not 0 < thermo.efficiency <= 1:
            raise ValueError(f"КПД ступени должен быть в диапазоне (0; 1], получено {thermo.efficiency}")

        L_ad = thermo.L_stage / thermo.efficiency
        c_ad = math.sqrt(2 * L_ad)
        a_crit_in = velocity_critical(k=K_GAS, R=R_GAS, T=thermo.T_in)
        lambda_ad = c_ad / a_crit_in
        pi_lambda_ad = pi_lambda(lam=lambda_ad, k=K_GAS)
        p_2_static_ad = pi_lambda_ad * thermo.p_in
        L_ad_stator = (1 - kin.reaction) * L_ad
        L_ad_rotor = kin.reaction * L_ad

        return TurbineStageWork(
            L_ad=L_ad, L_ad_stator=L_ad_stator, L_ad_rotor=L_ad_rotor, L_stage=thermo.L_stage,
            c_ad=c_ad, lambda_ad=lambda_ad, p_2_static_ad=p_2_static_ad, y=kin.u_mid / c_ad
        )

    def calculate_stator_outlet(
        self, flow: TurbineStageFlowParams, thermo: TurbineStageThermoParams,
        kin: TurbineStageKinematicParams, geom: TurbineStageGeometryParams,
        work: TurbineStageWork
    ) -> TurbineStatorOutlet:
        cp_gas = K_GAS / (K_GAS - 1) * R_GAS
        h_0_stator = (geom.inlet.tip - geom.inlet.hub) / 2
        h_1_stator = h_0_stator + (geom.h_2_rotor - h_0_stator) / 2

        c_1t = math.sqrt(2 * work.L_ad_stator)
        c_1 = kin.phi_cooling * c_1t
        T_1 = thermo.T_in - kin.phi_cooling ** 2 * work.L_ad_stator / cp_gas

        a_crit_in = velocity_critical(k=K_GAS, R=R_GAS, T=thermo.T_in)
        lambda_1t = c_1t / a_crit_in
        pi_lambda_1t = pi_lambda(lam=lambda_1t, k=K_GAS)
        p_1_static = pi_lambda_1t * thermo.p_in
        density_1 = p_1_static / (R_GAS * T_1)

        alpha_1_rad = math.asin(flow.G_in / (math.pi * geom.inlet.mid * h_1_stator * c_1 * density_1))
        alpha_1_deg = math.degrees(alpha_1_rad)

        w_1 = math.sqrt(c_1 ** 2 + kin.u_mid ** 2 - 2 * c_1 * kin.u_mid * math.cos(alpha_1_rad))
        beta_1_rad = math.asin(c_1 * math.sin(alpha_1_rad) / w_1)
        beta_1_deg = math.degrees(beta_1_rad)

        lambda_1 = c_1 / a_crit_in
        q_lambda_1 = q_lambda(lam=lambda_1, k=K_GAS)
        G_1 = flow.G_in - (flow.G_out - flow.G_in) / 2

        F_1 = area_calc(G=G_1, T=thermo.T_in, s=S_GAS, p=thermo.p_in, q=q_lambda_1 * math.sin(alpha_1_rad))
        section_1 = calculate_section_diameters(D_ref=geom.inlet.mid, F=F_1, mode_name=geom.mode, MODES_D=MODES_D)

        K_reaction = ((geom.inlet.mid / (geom.inlet.mid - h_1_stator) * math.cos(alpha_1_rad)) ** 2 + math.sin(alpha_1_rad) ** 2)
        reaction_hub = 1 - (1 - kin.reaction) * K_reaction

        triangle = VelocityTriangle(
            c=c_1, c_a=c_1 * math.sin(alpha_1_rad), c_u=c_1 * math.cos(alpha_1_rad),
            w=w_1, alpha_deg=alpha_1_deg, beta_deg=beta_1_deg
        )

        return TurbineStatorOutlet(
            triangle=triangle, T_static=T_1, p_static=p_1_static, density=density_1,
            lambda_1=lambda_1, lambda_1t=lambda_1t, q_lambda_1=q_lambda_1,
            area=F_1, section=section_1, h_stator_in=h_0_stator, h_stator_out=h_1_stator,
            K_reaction=K_reaction, reaction_hub=reaction_hub
        )

    def calculate_rotor_outlet(
        self, flow: TurbineStageFlowParams, thermo: TurbineStageThermoParams,
        kin: TurbineStageKinematicParams, geom: TurbineStageGeometryParams,
        work: TurbineStageWork, stator: TurbineStatorOutlet
    ) -> TurbineRotorOutlet:
        cp_gas = K_GAS / (K_GAS - 1) * R_GAS
        w_1 = stator.triangle.w
        c_1 = stator.triangle.c

        w_2 = kin.psi_cooling * math.sqrt(w_1 ** 2 + 2 * work.L_ad_rotor)
        T_w1 = thermo.T_in - (c_1 ** 2 - w_1 ** 2) / (2 * cp_gas)

        a_w1crit = velocity_critical(k=K_GAS, R=R_GAS, T=T_w1)
        lambda_w2 = w_2 / a_w1crit
        pi_lambda_w2 = pi_lambda(lam=lambda_w2, k=K_GAS)
        q_lambda_w2 = q_lambda(lam=lambda_w2, k=K_GAS)
        p_w2 = work.p_2_static_ad / pi_lambda_w2

        beta_2_rad = math.asin(flow.G_out * math.sqrt(T_w1) / (math.pi * geom.inlet.mid * geom.h_2_rotor * q_lambda_w2 * S_GAS * p_w2))
        beta_2_deg = math.degrees(beta_2_rad)

        c_2 = math.sqrt(w_2 ** 2 + kin.u_mid ** 2 - 2 * w_2 * kin.u_mid * math.cos(beta_2_rad))
        alpha_2_rad = math.asin(w_2 * math.sin(beta_2_rad) / c_2)
        alpha_2_deg = math.degrees(alpha_2_rad)

        T_2 = T_2 = thermo.T_in - work.L_stage / cp_gas
        a_2crit = velocity_critical(k=K_GAS, R=R_GAS, T=T_2)
        lambda_c2 = c_2 / a_2crit
        pi_lambda_c2 = pi_lambda(lam=lambda_c2, k=K_GAS)
        p_2 = work.p_2_static_ad / pi_lambda_c2
        q_lambda_2 = q_lambda(lam=lambda_c2, k=K_GAS)

        F_2 = area_calc(G=flow.G_out, T=T_2, s=S_GAS, p=p_2, q=q_lambda_2 * math.sin(alpha_2_rad))
        section_2 = calculate_section_diameters(D_ref=geom.inlet.mid, F=F_2, mode_name=geom.mode, MODES_D=MODES_D)

        triangle = VelocityTriangle(
            c=c_2, c_a=c_2 * math.sin(alpha_2_rad), c_u=c_2 * math.cos(alpha_2_rad),
            w=w_2, alpha_deg=alpha_2_deg, beta_deg=beta_2_deg
        )

        return TurbineRotorOutlet(
            triangle=triangle, T_relative_in=T_w1, T_total=T_2, p_relative_out=p_w2,
            p_total=p_2, lambda_w2=lambda_w2, lambda_c2=lambda_c2,
            q_lambda_w2=q_lambda_w2, q_lambda_2=q_lambda_2, area=F_2, section=section_2
        )

    def calculate_blade_rows(
        self, geom: TurbineStageGeometryParams, rows: TurbineStageBladeRowsParams,
        kin: TurbineStageKinematicParams, stator: TurbineStatorOutlet,
        rotor: TurbineRotorOutlet
    ) -> TurbineStageBladeRows:
        alpha_1_rad = math.radians(stator.triangle.alpha_deg)
        beta_1_rad = math.radians(stator.triangle.beta_deg)
        beta_2_rad = math.radians(rotor.triangle.beta_deg)

        delta_alpha = 180 - (90 + stator.triangle.alpha_deg)
        delta_beta = 180 - (stator.triangle.beta_deg + rotor.triangle.beta_deg)

        K_stator = 1 / math.sin(alpha_1_rad)
        K_rotor = math.sin(beta_1_rad) / math.sin(beta_2_rad)

        stator_grid = self.blade_row_grid_params(
            row=rows.stator, D_mean=geom.inlet.mid, K=K_stator, delta_angle_deg=delta_alpha,
            lambda_value=stator.lambda_1t, gamma_angle_deg=90 - stator.triangle.alpha_deg
        )

        rotor_grid = self.blade_row_grid_params(
            row=rows.rotor, D_mean=geom.inlet.mid, K=K_rotor, delta_angle_deg=delta_beta,
            lambda_value=rotor.lambda_w2 / kin.psi_cooling,
            gamma_angle_deg=stator.triangle.beta_deg - rotor.triangle.beta_deg
        )

        solidity_hub_rotor = math.pi * rotor.section.hub / rotor_grid.blade_count

        return TurbineStageBladeRows(
            stator_grid=stator_grid, rotor_grid=rotor_grid,
            delta_alpha_deg=delta_alpha, delta_beta_deg=delta_beta,
            K_stator=K_stator, K_rotor=K_rotor, solidity_hub_rotor=solidity_hub_rotor
        )

    def blade_row_grid_params(
        self, row: BladeRowInput, D_mean: float, K: float,
        delta_angle_deg: float, lambda_value: float, gamma_angle_deg: float
    ) -> BladeRowGridParams:
        def optimal_grid_step(K: float, delta_deg: float) -> float:
            delta_rad = math.radians(delta_deg)
            if K >= 1.5:
                return 0.327 / (K ** 0.327 * math.cbrt(delta_rad)) - 0.994 / K ** 0.385 + 1.314
            if K > 1:
                return (1.727 / K - 0.869) / math.cbrt(delta_rad) - 1.71 / K + 1.604
            raise ValueError(f"K must be > 1, got K = {K}")

        def optimal_grid_step_delta(lam: float) -> float:
            return -0.625 * lam ** 2 + 0.48 * lam + 0.016

        def gamma_calc(angle_deg: float) -> float:
            return 68.7 + 9.33e-4 * angle_deg - 6.052e-3 * angle_deg ** 2

        def trailing_edge_thickness_coeff(S_rel: float, t_rel: float) -> float:
            return 1 - 15 * S_rel ** 2 + (3.75 * t_rel - 0.6) * S_rel

        t_opt_rel = optimal_grid_step(K=K, delta_deg=delta_angle_deg)
        delta_t_opt_rel = optimal_grid_step_delta(lam=lambda_value)
        gamma_deg = gamma_calc(angle_deg=gamma_angle_deg)

        blade_chord = row.axial_chord / math.sin(math.radians(gamma_deg))
        trailing_edge_radius = row.trailing_edge_radius_rel * blade_chord
        S_out_rel = 2 * trailing_edge_radius / blade_chord
        K_crit = trailing_edge_thickness_coeff(S_rel=S_out_rel, t_rel=t_opt_rel)

        pitch = (1 + delta_t_opt_rel) * K_crit * t_opt_rel * blade_chord
        pitch_with_cooling = pitch * row.cooling_pitch_factor
        blade_count = round(math.pi * D_mean / pitch_with_cooling)
        solidity = blade_chord / pitch_with_cooling

        return BladeRowGridParams(
            t_opt_rel=t_opt_rel, delta_t_opt_rel=delta_t_opt_rel, gamma_deg=gamma_deg,
            blade_chord=blade_chord, trailing_edge_radius=trailing_edge_radius, S_out_rel=S_out_rel,
            K_crit=K_crit, pitch=pitch, pitch_with_cooling=pitch_with_cooling,
            solidity=solidity, blade_count=blade_count
        )

    def assemble_result(
        self, thermo: TurbineStageThermoParams, work: TurbineStageWork,
        stator: TurbineStatorOutlet, rotor: TurbineRotorOutlet,
        blade_rows: TurbineStageBladeRows
    ) -> TurbineStageResult:
        geom = self.params.geometry
        kin = self.params.kinematics

        thermodynamics = TurbineStageThermodynamics(
            L_ad=work.L_ad, L_ad_stator=work.L_ad_stator, L_ad_rotor=work.L_ad_rotor,
            L_stage=work.L_stage, T_in=thermo.T_in, T_1=stator.T_static,
            T_w1=rotor.T_relative_in, T_2=rotor.T_total, T_out_target=thermo.T_out,
            p_in=thermo.p_in, p_1_static=stator.p_static, p_2_static=work.p_2_static_ad,
            p_w2=rotor.p_relative_out, p_out=rotor.p_total, pi_stage=thermo.pi
        )

        geometry = TurbineStageGeometry(
            inlet=geom.inlet, stator_outlet=stator.section, rotor_outlet=rotor.section,
            h_0_stator=stator.h_stator_in, h_1_stator=stator.h_stator_out,
            h_2_rotor=geom.h_2_rotor, F_1=stator.area, F_2=rotor.area,
        )

        loading = TurbineStageLoading(
            reaction_mean=kin.reaction, reaction_hub=stator.reaction_hub,
            K_reaction=stator.K_reaction, delta_alpha_deg=blade_rows.delta_alpha_deg,
            delta_beta_deg=blade_rows.delta_beta_deg, K_stator=blade_rows.K_stator,
            K_rotor=blade_rows.K_rotor
        )

        velocity = TurbineStageVelocity(
            stator_outlet=stator.triangle, rotor_outlet=rotor.triangle,
            c_ad=work.c_ad, lambda_ad=work.lambda_ad, y=work.y,
            lambda_1=stator.lambda_1, lambda_1t=stator.lambda_1t,
            lambda_w2=rotor.lambda_w2, lambda_c2=rotor.lambda_c2
        )

        return TurbineStageResult(
            thermodynamics=thermodynamics, geometry=geometry, velocity=velocity,
            loading=loading, stator_grid=blade_rows.stator_grid,
            rotor_grid=blade_rows.rotor_grid,
            solidity_hub_rotor=blade_rows.solidity_hub_rotor
        )
