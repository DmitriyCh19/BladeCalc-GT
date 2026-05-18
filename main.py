from machine.hpt import HPTParameters, HighPressureTurbine
from machine.hpc import HPCParameters, HighPressureCompressor
from visualization.flowpath_builders import build_hpc_flowpath_data, build_hpt_flowpath_data
from visualization.flowpath_plot import plot_machine_flowpath
from core.gas_func import velocity_critical
from configs.constants import *

hpt_params = HPTParameters(
    mode='mid',

    T_gas=1640,
    p_gas=2877336,

    T_out=1208.56,
    p_out=597099,

    T_cooling_air=862.42,

    z=2,
    efficiency=0.87,
    y=0.55,
    theta=0.45,

    G_in=80.3982,
    G_cooling_rel=0.06,
    G_hpc_in=85.4106,

    density=8100,
    blade_force_coeff=0.7,
    sigma_allow=380,
    k_sigma=1.8,

    alpha_2_deg=65,

    lambda_2=0.45,
    lambda_gas=0.25,

    L_stage=542448,
)

hpc_params = HPCParameters(
    mode='tip',

    T_in=378.76,
    T_out=862.42,
    p_out=3060995,
    pi=13.6,
    d_hub_in_rel=0.74,
    K_g=0.96,
    efficiency=0.816,
    lambda_in=0.6,
    lambda_out=0.3,

    G_in=83.539,
    G_cooling_rel=0.06
)


if __name__ == '__main__':
    K_gg = 0.45

    hpt = HighPressureTurbine(hpt_params)
    hpc = HighPressureCompressor(hpc_params)
    a_crit_in = velocity_critical(k=K_AIR, R=R_AIR, T=hpc_params.T_in)
    c_in = a_crit_in * hpc_params.lambda_in
    print(c_in)

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
    print(f'u внешний диаметр КВД 1 ступень ≤ 550 м/с = {hpc.u_tip_in}')
    print(f'частота вращения РВД {hpt.n_hp}')
    print(z_hpc)

    # stages_result = hpc.calculate_stages(
    #     z=9,
    #     L_stage_rel=[0.24227, 0.252, 0.262, 0.271, 0.281, 0.271, 0.261, 0.251, 0.241],
    #     eff_stage=[0.87, 0.89, 0.89, 0.9, 0.91, 0.91, 0.9, 0.885, 0.885],
    #     eff_rel=[1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.935, 0.93],
    #     reaction_stage=[0.5, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57],
    #     c1a_stage=[193.4903, 191.9, 190.3, 188.7, 187.2, 187.2, 185.6, 184.01, 182.4, 180.8988],
    #     h_rot_rel=[3.5, 3.5, 3.3, 3.1, 2.9, 2.7, 2.4, 2.2, 2],
    #     eff_rotor=0.93,
    # )

    # for i, stage in enumerate(hpc.stage_results, start=1):
    #     print(f'\nСтупень {i}')
    #     print(f'T_in = {stage.thermodynamics.T_in:.2f} К')
    #     print(f'T_out = {stage.thermodynamics.T_out:.2f} К')
    #     print(f'p_in = {stage.thermodynamics.p_in:.0f} Па')
    #     print(f'p_out = {stage.thermodynamics.p_out:.0f} Па')
    #     print(f'pi_stage = {stage.thermodynamics.pi:.4f}')
    #     print(f'Ротор = {stage.rotor}')
    #     print(f'НА = {stage.stator}')

    # print('\nИтог по компрессору')
    # print(f'T_out = {stages_result.T_out:.2f} К')
    # print(f'p_out = {stages_result.p_out:.0f} Па')
    # print(f'pi_total = {stages_result.pi_total:.4f}')
    # print(f'L_total = {stages_result.L_total:.2f} Дж/кг')

    # print(
    # f'Ступень {i}: РК = {stage.length.rotor:.4f} м, '
    # f'НА = {stage.length.stator:.4f} м, '
    # f'зазор = {stage.length.gap:.4f} м, '
    # f'всего = {stage.length.total:.4f} м'
    # )

    # print('========================================')


    # hpt_stages = hpt.calculate_stages(
    # L_stage=[393673],
    # eff_stage=[0.88],
    # reaction_stage=[0.32],
    # stator_axial_chord=[0.0478],
    # rotor_axial_chord=[0.0427],
    # stator_trailing_edge_radius_rel=[0.015],
    # rotor_trailing_edge_radius_rel=[0.015],
    # cooling_rel_stage=[0.06],
    # phi_stage=[0.965],
    # psi_stage=[0.95],
    # stator_cooling_pitch_factor=[0.95],
    # rotor_cooling_pitch_factor=[1.1],
    # )

    # print('\nИтог по ТВД')
    # print(f'T_out = {hpt_stages.T_out:.2f} К')
    # print(f'p_out = {hpt_stages.p_out:.0f} Па')
    # print(f'pi_total = {hpt_stages.pi_total:.4f}')
    # print(f'L_total = {hpt_stages.L_total:.2f} Дж/кг')

    # for i, stage in enumerate(hpt_stages.stages, start=1):
    #     print(f'\nСтупень турбины {i}')
    #     print(f'T_in / T_out = {stage.thermodynamics.T_in:.2f} / {stage.thermodynamics.T_2:.2f} К')
    #     print(f'p_in / p_out = {stage.thermodynamics.p_in:.0f} / {stage.thermodynamics.p_out:.0f} Па')
    #     print(f'α1 / β1 = {stage.velocity.stator_outlet.alpha_deg:.2f} / {stage.velocity.stator_outlet.beta_deg:.2f} град')
    #     print(f'α2 / β2 = {stage.velocity.rotor_outlet.alpha_deg:.2f} / {stage.velocity.rotor_outlet.beta_deg:.2f} град')
    #     print(f'СА: z = {stage.stator_grid.blade_count}, РК: z = {stage.rotor_grid.blade_count}')
    #     print(f'y = {stage.velocity.y:.4f}')

    # print('\nДлина ТВД')
    # print(f'L_total = {hpt_stages.length_total:.4f} м')

    # for i, stage in enumerate(hpt_stages.stages, start=1):
    #     print(
    #         f'Ступень {i}: СА = {stage.length.stator:.4f} м, '
    #         f'РК = {stage.length.rotor:.4f} м, '
    #         f'зазор = {stage.length.gap:.4f} м, '
    #         f'всего = {stage.length.total:.4f} м'
    #     )

    # print('\nКонтроль совпадения выхода КВД и последней ступени')
    # for name, value in stages_result.checks.items():
    #     print(f'{name}: {value:.4f} %')

    # print('\nКонтроль совпадения выхода ТВД и последней ступени')
    # for name, value in hpt_stages.checks.items():
    #     print(f'{name}: {value:.4f} %')




    # hpc_plot_data = build_hpc_flowpath_data(hpc)
    # plot_machine_flowpath(hpc_plot_data)
    # hpt_plot_data = build_hpt_flowpath_data(hpt)
    # plot_machine_flowpath(hpt_plot_data)