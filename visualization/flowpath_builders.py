from visualization.flowpath_models import (
    FlowPathSection,
    FlowPathRow,
    FlowPathStation,
    MachineFlowPathPlotData,
)


def build_hpc_flowpath_data(hpc) -> MachineFlowPathPlotData:
    if hpc.geometry is None:
        raise ValueError("Перед построением проточной части нужно вызвать hpc.calculate(n).")
    if not hpc.stage_results:
        raise ValueError("Перед построением проточной части нужно рассчитать ступени КВД.")

    sections = []
    rows = []
    stations = []
    x = 0.0

    first_stage = hpc.stage_results[0]
    
    sections.append(
        FlowPathSection(
            x=x,
            hub=first_stage.rotor.inlet.hub / 2,
            tip=first_stage.rotor.inlet.tip / 2,
        )
    )
    
    stations.append(
        FlowPathStation(
            name="Вход КВД",
            x=x,
            T_total=first_stage.thermodynamics.T_in,
            p_total=first_stage.thermodynamics.p_in,
        )
    )

    for i, stage in enumerate(hpc.stage_results, start=1):
        length = stage.length

        if length is None:
            raise ValueError(f"Для ступени КВД {i} не рассчитана длина.")

        # Зазор перед РК: канал идет к входу РК
        x += length.gap
        sections.append(
            FlowPathSection(
                x=x,
                hub=stage.rotor.inlet.hub / 2,
                tip=stage.rotor.inlet.tip / 2,
            )
        )

        # РК компрессора: вход РК -> выход РК
        rows.append(
            FlowPathRow(
                machine="HPC",
                stage_index=i,
                row_type="rotor",
                x0=x,
                x1=x + length.rotor,
                hub_in=stage.rotor.inlet.hub / 2,
                hub_out=stage.rotor.outlet.hub / 2,
                tip_in=(stage.rotor.inlet.hub / 2 + stage.rotor.blade_height_in * 0.95),
                tip_out=(stage.rotor.outlet.hub / 2 + stage.rotor.blade_height_out * 0.95),
            )
        )

        x += length.rotor
        sections.append(
            FlowPathSection(
                x=x,
                hub=stage.rotor.outlet.hub / 2,
                tip=stage.rotor.outlet.tip / 2,
            )
        )

        stations.append(
            FlowPathStation(
                name=f"КВД ст. {i} после РК",
                x=x,
                T_total=stage.thermodynamics.T_out,
                p_total=stage.thermodynamics.p_2,
            )
        )

        # Зазор между РК и НА: геометрия на входе НА та же, что на выходе РК
        x += length.gap
        sections.append(
            FlowPathSection(
                x=x,
                hub=stage.rotor.outlet.hub / 2,
                tip=stage.rotor.outlet.tip / 2,
            )
        )

        if i < len(hpc.stage_results):
            next_stage = hpc.stage_results[i]
            stator_out = next_stage.rotor.inlet
        else:
            stator_out = hpc.geometry.outlet

        # НА компрессора: D_in = D_out текущего РК, D_out = D_in следующего РК
        rows.append(
            FlowPathRow(
                machine="HPC",
                stage_index=i,
                row_type="stator",
                x0=x,
                x1=x + length.stator,
                hub_in=stage.rotor.outlet.hub / 2,
                hub_out=stator_out.hub / 2,
                tip_in=stage.rotor.outlet.tip / 2,
                tip_out=stator_out.tip / 2,
            )
        )

        x += length.stator
        sections.append(
            FlowPathSection(
                x=x,
                hub=stator_out.hub / 2,
                tip=stator_out.tip / 2,
            )
        )

        # В НА полные параметры не меняются
        stations.append(
            FlowPathStation(
                name=f"КВД ст. {i} после НА",
                x=x,
                T_total=stage.thermodynamics.T_out,
                p_total=stage.thermodynamics.p_2,
            )
        )

    return MachineFlowPathPlotData(
        name="КВД",
        sections=sections,
        rows=rows,
        stations=stations,
        x_max=x,
    )


def build_hpt_flowpath_data(hpt) -> MachineFlowPathPlotData:
    if hpt.geometry is None:
        raise ValueError("Перед построением проточной части нужно вызвать hpt.calculate().")
    if not hpt.stage_results:
        raise ValueError("Перед построением проточной части нужно рассчитать ступени ТВД.")

    sections = []
    rows = []
    stations = []
    x = 0.0

    first_stage = hpt.stage_results[0]

    sections.append(
        FlowPathSection(
            x=x,
            hub=hpt.geometry.inlet.hub / 2,
            tip=hpt.geometry.inlet.tip / 2,
        )
    )

    stations.append(
        FlowPathStation(
            name="Вход ТВД",
            x=x,
            T_total=first_stage.thermodynamics.T_in,
            p_total=first_stage.thermodynamics.p_in,
        )
    )

    for i, stage in enumerate(hpt.stage_results, start=1):
        length = stage.length

        if length is None:
            raise ValueError(f"Для ступени ТВД {i} не рассчитана длина.")

        # СА турбины: вход ступени -> выход СА
        rows.append(
            FlowPathRow(
                machine="HPT",
                stage_index=i,
                row_type="stator",
                x0=x,
                x1=x + length.stator,
                hub_in=stage.geometry.inlet.hub / 2,
                hub_out=stage.geometry.stator_outlet.hub / 2,
                tip_in=stage.geometry.inlet.tip / 2,
                tip_out=stage.geometry.stator_outlet.tip / 2,
            )
        )

        x += length.stator
        sections.append(
            FlowPathSection(
                x=x,
                hub=stage.geometry.stator_outlet.hub / 2,
                tip=stage.geometry.stator_outlet.tip / 2,
            )
        )

        # В СА полные параметры не меняются
        stations.append(
            FlowPathStation(
                name=f"ТВД ст. {i} после СА",
                x=x,
                T_total=stage.thermodynamics.T_in,
                p_total=stage.thermodynamics.p_in,
            )
        )

        # Зазор между СА и РК
        x += length.gap
        sections.append(
            FlowPathSection(
                x=x,
                hub=stage.geometry.stator_outlet.hub / 2,
                tip=stage.geometry.stator_outlet.tip / 2,
            )
        )

        # РК турбины: выход СА / вход РК -> выход РК
        rows.append(
            FlowPathRow(
                machine="HPT",
                stage_index=i,
                row_type="rotor",
                x0=x,
                x1=x + length.rotor,
                hub_in=stage.geometry.stator_outlet.hub / 2,
                hub_out=stage.geometry.rotor_outlet.hub / 2,
                tip_in=(stage.geometry.stator_outlet.hub / 2 +  0.98 *(stage.geometry.stator_outlet.tip - stage.geometry.stator_outlet.hub) / 2),
                tip_out=(stage.geometry.rotor_outlet.hub / 2 +  0.98 *(stage.geometry.rotor_outlet.tip - stage.geometry.rotor_outlet.hub) / 2),
            )
        )

        x += length.rotor
        sections.append(
            FlowPathSection(
                x=x,
                hub=stage.geometry.rotor_outlet.hub / 2,
                tip=stage.geometry.rotor_outlet.tip / 2,
            )
        )

        stations.append(
            FlowPathStation(
                name=f"ТВД ст. {i} после РК",
                x=x,
                T_total=stage.thermodynamics.T_2,
                p_total=stage.thermodynamics.p_out,
            )
        )

        # Зазор после РК
        x += length.gap
        sections.append(
            FlowPathSection(
                x=x,
                hub=stage.geometry.rotor_outlet.hub / 2,
                tip=stage.geometry.rotor_outlet.tip / 2,
            )
        )

    return MachineFlowPathPlotData(
        name="ТВД",
        sections=sections,
        rows=rows,
        stations=stations,
        x_max=x,
    )