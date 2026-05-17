from machine.hpt import HPTParameters, HighPressureTurbine
from machine.hpc import HPCParameters, HighPressureCompressor
from machine.turbine_stage import TurbineStage
from machine.turbine_stage_models import (TurbineStageParameters, BladeRowInput, TurbineStageBladeRowsParams, 
TurbineStageFlowParams, TurbineStageGeometryParams, TurbineStageKinematicParams, TurbineStageThermoParams)
from core.geometry_models import SectionDiameters


z = 9
alpha11_deg = 60.546
L_stage_rel = [0.24227, 0.252, 0.262, 0.271, 0.281, 0.271, 0.261, 0.251, 0.241]
D_tip = 0.5923
eff_stage = [0.87, 0.89, 0.89, 0.9, 0.91, 0.91, 0.9, 0.885, 0.885]
eff_rel = [1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.935, 0.93]
reaction_stage = [0.5, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57]
c1a_stage = [193.4903, 191.9, 190.3, 188.7, 187.2, 187.2, 185.6, 184.01, 182.4]
h_rot_rel = [3.5, 3.5, 3.3, 3.1, 2.9, 2.7, 2.4, 2.2, 2]
G = 72.4359
d_hub_rel = 0.74
lambda_in = 0.58
u_tip_1 = 393.4265
n = 12692
T_in = 437.9788
p_in = 355623
pi = 6.58
D_1mid_1 = 0.521
D_1hub_1 = 0.4383

hpt_params = HPTParameters(
    mode='mid',

    T_gas=1665,
    p_gas=2199900,

    T_out=1327.36,
    p_out=765440,

    T_cooling_air=796.99,

    z=1,
    efficiency=0.88,
    y=0.445,
    theta=0.58,

    G_in=66.496,
    G_cooling_rel=0.06,
    G_hpc_in=72.435,

    density=8100,
    blade_force_coeff=0.7,
    sigma_allow=380,
    k_sigma=1.8,

    alpha_2_deg=65,

    lambda_2=0.497,
    lambda_gas=0.25,

    L_stage=393673,
)

hpc_params = HPCParameters(
    mode='tip',

    T_in=437.97,
    T_out=796.99,
    p_out=2.3402e6,
    pi=6.58,
    d_hub_in_rel=0.74,
    K_g=0.96,
    efficiency=0.87,
    lambda_in=0.58,
    lambda_out=0.35,

    G_in=72.435,
    G_cooling_rel=0.1
)

turbin_stage_par = TurbineStageParameters(
    flow=TurbineStageFlowParams(
        G_in=66.496,
        G_out=70.84,
    ),
    thermo=TurbineStageThermoParams(
        p_in=2199800,
        pi=2.8739,
        T_in=1665,
        T_out=1327.36,
        L_stage=393673,
    ),
    kinematics=TurbineStageKinematicParams(
        u_mid=420.922,
        y=0.445,
        reaction=0.32,
        phi_cooling=0.965,
        psi_cooling=0.95,
    ),
    geometry=TurbineStageGeometryParams(
        mode='mid',
        inlet=SectionDiameters(hub=0.5919, mid=0.6338, tip=0.6731),
        h_2_rotor=0.0666,
    ),
    blade_rows=TurbineStageBladeRowsParams(
        stator=BladeRowInput(
            axial_chord=0.0478,
            trailing_edge_radius_rel=0.015,
            cooling_pitch_factor=0.95,
        ),
        rotor=BladeRowInput(
            axial_chord=0.0427,
            trailing_edge_radius_rel=0.015,
            cooling_pitch_factor=1.1,
        ),
    ),
)

if __name__ == '__main__':
    K_gg = 0.4

    hpt = HighPressureTurbine(hpt_params)
    hpc = HighPressureCompressor(hpc_params)

    hpt.calculate()
    hpc.calculate(n=hpt.n_hp)
    z_hpc = (hpt.D_mean / hpc.D_mean) ** 2 * hpt.params.z / K_gg ** 2

    print(hpt.geometry)
    print(f'u = {hpt.u_mid} ')
    # print(f'Относительная высота лопатки = {hpt.blade_height_rel}')
    # print(f'Высота лопатки на выходе из ТВД ≥0.025м = {hpt.h_blade_1}')
    print(f'Высота лопатки на выходе из ТВД ≥0.025м = {hpt.h_blade_out}')
    print(hpc.geometry)
    print(f'c_out КВД 150...180 = {hpc.c_out}')
    print(f'u внешний диаметр КВД 1 ступень ≤ 550 м/с = {hpc.u_tip_in} ')
    print(z_hpc)

    c1a_stage.append(180.8988)

    stages_result = hpc.calculate_stages(
        z=len(L_stage_rel),
        alpha1_deg=alpha11_deg,
        L_stage_rel=L_stage_rel,
        eff_stage=eff_stage,
        eff_rel=eff_rel,
        reaction_stage=reaction_stage,
        c1a_stage=c1a_stage,
        h_rot_rel=h_rot_rel,

        D_const=D_tip,
        D_tip=D_tip,
        G=G,
        d_hub_rel=d_hub_rel,
        lambda_in=lambda_in,
        u_tip=u_tip_1,
        T_in=T_in,
        p_in=p_in,
        D_1mid=D_1mid_1,
        D_1hub=D_1hub_1,
        K_g=0.96,
        eff_rotor=0.93,
        mode='tip',
    )

    for i, stage in enumerate(hpc.stage_results, start=1):
        print(f'\nСтупень {i}')
        print(f'T_in = {stage.thermodynamics.T_in:.2f} К')
        print(f'T_out = {stage.thermodynamics.T_out:.2f} К')
        print(f'p_in = {stage.thermodynamics.p_in:.0f} Па')
        print(f'p_out = {stage.thermodynamics.p_out:.0f} Па')
        print(f'pi_stage = {stage.thermodynamics.pi:.4f}')
        print(f'Ротор = {stage.rotor}')
        print(f'НА = {stage.stator}')

    print('\nИтог по компрессору')
    print(f'T_out = {stages_result.T_out:.2f} К')
    print(f'p_out = {stages_result.p_out:.0f} Па')
    print(f'pi_total = {stages_result.pi_total:.4f}')
    print(f'L_total = {stages_result.L_total:.2f} Дж/кг')

    print('========================================')

    turb_stage = TurbineStage(params=turbin_stage_par)
    turb_result = turb_stage.calculate()

    print('\n' + '=' * 60)
    print('КОНТРОЛЬ РАСЧЁТА СТУПЕНИ ТУРБИНЫ')
    print('=' * 60)

    print(f'Lад = {turb_result.thermodynamics.L_ad:.2f} Дж/кг')
    print(f'T2 расчётная / заданная = {turb_result.thermodynamics.T_2:.2f} / {turb_result.thermodynamics.T_out_target:.2f} К')
    print(f'p2 расчётное = {turb_result.thermodynamics.p_out:.0f} Па')
    print(f'α1 / β1 = {turb_result.velocity.stator_outlet.alpha_deg:.2f} / {turb_result.velocity.stator_outlet.beta_deg:.2f} град')
    print(f'α2 / β2 = {turb_result.velocity.rotor_outlet.alpha_deg:.2f} / {turb_result.velocity.rotor_outlet.beta_deg:.2f} град')
    print(f'F1 / F2 = {turb_result.geometry.F_1:.5f} / {turb_result.geometry.F_2:.5f} м²')
    print(f'ρ ср / ρ втулки = {turb_result.loading.reaction_mean:.4f} / {turb_result.loading.reaction_hub:.4f}')
    print(f'СА: z = {turb_result.stator_grid.blade_count}, b = {turb_result.stator_grid.blade_chord:.5f} м, σ = {turb_result.stator_grid.solidity:.4f}')
    print(f'РК: z = {turb_result.rotor_grid.blade_count}, b = {turb_result.rotor_grid.blade_chord:.5f} м, σ = {turb_result.rotor_grid.solidity:.4f}')
    print('=' * 60)