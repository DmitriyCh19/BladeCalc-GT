import math
from dataclasses import dataclass

from configs.constants import *
from configs.modes import *
from core.gas_func import velocity_critical, q_lambda, local_speed_velocity
from core.geometry import area_calc, relative_diameter_hub, calculate_section_diameters
from core.geometry_models import VelocityTriangle, RotorGeometry, StatorGeometry

@dataclass
class CompressorStageParameters:
    mode: str
    D_const: float

    L_rel: float
    efficiency: float
    eff_tip_rel: float
    reaction: float
    c_1a: float
    c_3a: float
    h_rot_rel: float
    G: float
    u_tip: float
    T_in: float
    p_in: float
    K_g: float
    eff_rotor: float

@dataclass
class StageThermodynamics:
    L_stage: float

    T_in: float
    T_out: float

    p_in: float
    p_2:float
    p_out: float

    pi: float

@dataclass
class InletSolution:
    area: float

    u_mid: float

    lambda_1: float

    L_ku_rel: float

@dataclass
class CompressorStageResult:
    thermodynamics: StageThermodynamics
    inlet_triangle: VelocityTriangle
    outlet_triangle: VelocityTriangle
    rotor: RotorGeometry

    u_mid: float
    flow_coefficient: float
    J: float
    lambda_1: float
    lambda_2: float
    mach_1_rel: float
    mach_2_abs: float
    F_inlet: float
    F_outlet: float
    a_inlet: float
    a_outlet: float

    stator: StatorGeometry | None = None

class CompressorStage:
    def __init__(self, params:CompressorStageParameters):
        self.params = params

        self.inlet_triangle = None
        self.outlet_triangle = None

        self.rotor = None
        self.stage = None

    def calculate(self):
        thermo = self.calculate_thermodynamics()
        self.calculate_rotor(thermo=thermo)

    def calculate_rotor(self, thermo:StageThermodynamics):    
        L_ku = thermo.L_stage* self.params.eff_tip_rel

        # u_mid_in = self.params.u_tip * math.sqrt((1 + self.params.d_hub_rel ** 2) / 2)
        a_1crit = velocity_critical(k=K_AIR, R=R_AIR, T=thermo.T_in)

        # пропущенны пункты 9-10

        lambda_1a = self.params.c_1a / a_1crit
        q_lambda_1a = q_lambda(lam=lambda_1a, k=K_AIR)
        
        inlet = self.solve_inlet_geometry(q_lambda_1a, L_ku, a_1crit)
        u_mid = inlet.u_mid
        flow_coefficient = self.params.c_1a / u_mid

        stage_1 = calculate_section_diameters(
            D_ref=self.params.D_const,
            F=inlet.area,
            mode_name=self.params.mode,
            MODES_D=MODES_D
        )

        h_1 = (stage_1.tip - stage_1.hub) / 2

        alpha1_rad = math.radians(self.inlet_triangle.alpha_deg)
        self.inlet_triangle.beta_deg = math.degrees(math.atan(flow_coefficient / (1 - flow_coefficient * (1 / math.tan(alpha1_rad)))))

        self.inlet_triangle.w = math.sqrt(self.inlet_triangle.c_a ** 2 + (u_mid - self.inlet_triangle.c_u) ** 2)

        a_1 = local_speed_velocity(a_crit=a_1crit, k=K_AIR, lam=inlet.lambda_1)

        mach_1 = self.inlet_triangle.w  / a_1

        par_L_c = inlet.L_ku_rel / flow_coefficient
        par_rho_c = self.params.reaction / flow_coefficient
        par_L_c_bt = 0.7 - 0.27 * par_rho_c + 0.16 * par_rho_c ** 2
        J = par_L_c / par_L_c_bt
        rotor_solidity = 0.225 + 0.275 * J + 0.5 * J ** 2

        z_rotor = self.params.h_rot_rel * rotor_solidity * math.pi * stage_1.mid / h_1
        z_rotor = round(z_rotor)
        h_rot_rel = h_1 * z_rotor / (rotor_solidity * math.pi * stage_1.mid)

        rotor_blade_chord = h_1 / h_rot_rel
        
        c_2u = u_mid * ((1 - self.params.reaction) + inlet.L_ku_rel / 2)
        c_2a = (self.params.c_1a + self.params.c_3a) / 2

        self.outlet_triangle = VelocityTriangle(c=math.sqrt(c_2a ** 2 + c_2u ** 2),
                                                c_a=c_2a,
                                                c_u=c_2u)
        
        a_2crit = velocity_critical(k=K_AIR, R=R_AIR, T=thermo.T_out)
        lambda_2 = self.outlet_triangle.c / a_2crit
        a_2 = local_speed_velocity(a_crit=a_2crit, k=K_AIR, lam=lambda_2)

        alpha2_rad = math.asin(self.outlet_triangle.c_a / self.outlet_triangle.c)
        mach_2 = self.outlet_triangle.c / a_2
        self.outlet_triangle.alpha_deg = math.degrees(alpha2_rad)
        
        sigma_stator = thermo.p_out / thermo.p_2

        q_lambda_2 = q_lambda(lam=lambda_2, k=K_AIR)
        F_outlet = area_calc(G=self.params.G, T=thermo.T_out, s=S_AIR, p=thermo.p_2,
                            q=(q_lambda_2 * self.params.K_g * math.sin(alpha2_rad)))
        
        d_2_hub_rel = relative_diameter_hub(D_ref=self.params.D_const, F=F_outlet, mode_name=self.params.mode)

        stage_2 = calculate_section_diameters(
            D_ref=self.params.D_const,
            F=F_outlet,
            mode_name=self.params.mode,
            MODES_D=MODES_D
        )
        h_2 = (stage_2.tip - stage_2.hub) / 2
        self.outlet_triangle.w = math.sqrt(self.outlet_triangle.c_a ** 2 + (u_mid - self.outlet_triangle.c_u) ** 2)

        beta_2_rad = math.asin(self.outlet_triangle.c_a / self.outlet_triangle.w)
        self.outlet_triangle.beta_deg = math.degrees(beta_2_rad)

        delta_beta = self.outlet_triangle.beta_deg - self.inlet_triangle.beta_deg

        self.rotor = RotorGeometry(inlet=stage_1, 
                              outlet=stage_2, 
                              blade_height_in=h_1,
                              blade_height_out=h_2,
                              blade_chord=rotor_blade_chord,
                              blade_count=z_rotor, solidity=rotor_solidity)
        
        self.stage = CompressorStageResult(thermodynamics=thermo,
                                           inlet_triangle=self.inlet_triangle,
                                           outlet_triangle=self.outlet_triangle,
                                           rotor=self.rotor,
                                           u_mid=u_mid,
                                           flow_coefficient=flow_coefficient,
                                           J=J, lambda_1=inlet.lambda_1, lambda_2=lambda_2,
                                           mach_1_rel=mach_1, mach_2_abs=mach_2,
                                           F_inlet=inlet.area, F_outlet=F_outlet,
                                           a_inlet=a_1, a_outlet=a_2)

    def print_stage_result(self) -> None:
        # print(f'u_mid = {u_mid}')
        # print(f'a_crit = {a_1crit}')
        # print(f'lambda 1a = {lambda_1a}')
        # print(f'q_lambda_a1 {q_lambda_1a}')
        # print(f'u_mid = {u_mid}')
        # print(f'коэфф расхода {flow_coefficient}')
        # print(f'Диаметры на входе {stage_1}')
        # print(f'высота лопатки РК {h_1}')
        # print(f'скорость звука на входе {a_1}')
        print(f'Треугольник скоростей вход {self.inlet_triangle}')
        # print(f'безразмерный параметр J {J}')
        # print(f'густота решетки РК {rotor_solidity}')
        # print(f'число лопаток РК {z_rotor}')
        # print(f'относительная высота рабочей лопатки {h_rot_rel}')
        # print(f'хорда рабочих лопаток {rotor_blade_chord}')
        # print(f'lamda 2 = {lambda_2}')
        # print(f'скорость звука на выходе из РК {a_2}')
        print(f'Треугольник выхода из РК {self.outlet_triangle}')
        print(f'Геометрия ротора {self.rotor}')
        # print(f'Площадь на выходе из РК {F_outlet}')
        # print(f'относительный диаметр втулки на выходе из РК {d_2_hub_rel}')
        # print(f'Диаметры на выходе из РК {stage_2}')
        # print(f'высота лопатки на выходе из РК {h_2}')

    def calculate_thermodynamics(self) -> StageThermodynamics:
        L_st = self.params.L_rel * self.params.u_tip ** 2
        T_out = self.params.T_in + L_st / (K_AIR / (K_AIR - 1) * R_AIR) # T_in - температура перед прошлой ступенью
        pi = (L_st * self.params.efficiency / (K_AIR / (K_AIR - 1) * R_AIR * self.params.T_in) + 1) ** (K_AIR / (K_AIR - 1))
        p_out = self.params.p_in * pi
        p_2 = self.params.p_in * (1 + L_st * self.params.eff_rotor / (K_AIR / (K_AIR - 1) * R_AIR * self.params.T_in)) ** (K_AIR / (K_AIR - 1))
        return StageThermodynamics(L_stage=L_st, T_in=self.params.T_in,
                                   T_out=T_out, p_in=self.params.p_in, p_2=p_2, p_out=p_out, pi=pi)
    
    def solve_inlet_geometry(self, q_lam_1a, L_ku, a_crit) -> InletSolution: 
        delta = 10
        alpha_deg = 90
        q_lambda_1 = q_lam_1a
        while delta > 1:
            alpha_old = alpha_deg
            alpha_rad = math.radians(alpha_old)
            F = area_calc(G=self.params.G, T=self.params.T_in, s=S_AIR, p=self.params.p_in,
                            q=(q_lambda_1 * self.params.K_g * math.sin(alpha_rad)))

            d_hub_rel = relative_diameter_hub(D_ref=self.params.D_const, F=F, mode_name=self.params.mode)

            u_mid = self.params.u_tip * math.sqrt((1 + d_hub_rel ** 2) / 2) # реализовать расчет при разных D=const

            L_ku_rel = L_ku / u_mid ** 2

            c_u = u_mid * ((1 - self.params.reaction) - L_ku_rel / 2)

            c = math.sqrt(self.params.c_1a ** 2 + c_u ** 2)

            lambda_1 = c / a_crit
            q_lambda_1 = q_lambda(lam=lambda_1, k=K_AIR)

            alpha_deg = math.degrees(math.asin(self.params.c_1a / c))
            delta = abs((alpha_deg - alpha_old) / alpha_deg) * 100
            
        self.inlet_triangle = VelocityTriangle(
            c=c,
            c_a=self.params.c_1a,
            c_u=c_u,
            alpha_deg=alpha_deg
        )
        return InletSolution(area=F, u_mid=u_mid, lambda_1=lambda_1, L_ku_rel=L_ku_rel)
    
    def calculate_stator(self, alpha_out_deg: float):
        delta_alpha = alpha_out_deg - self.outlet_triangle.alpha_deg

        alpha_rel = (alpha_out_deg / 100)
        delta_alpha_b_t = 0.037 + 0.1 * alpha_rel + 0.262 * alpha_rel ** 2

        E = delta_alpha / 100 / delta_alpha_b_t
        if E <= 1:
            stator_solidity = 0.231 - 0.135 * E + 0.909 * E ** 2
        else:
            stator_solidity = 10 * (0.981 - 1.788 * E + 0.912 * E ** 2)
        z_stator = self.params.h_rot_rel * stator_solidity * math.pi * self.rotor.outlet.mid / self.rotor.blade_height_out
        z_stator = round(z_stator)
        h_stat_rel = z_stator * self.rotor.blade_height_out / (stator_solidity * math.pi * self.rotor.outlet.mid)

        stator_blade_chord = self.rotor.blade_height_out / h_stat_rel
        self.stator = StatorGeometry(blade_chord=stator_blade_chord, blade_count=z_stator, solidity=stator_solidity)

        self.stage.stator = self.stator



if __name__ == '__main__':
    pass