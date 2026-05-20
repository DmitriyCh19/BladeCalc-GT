import math


def distribute_efficiency_normal(
        eta_total: float,
        stage_count: int,
        *,
        spread: float = 0.02,
        sigma_rel: float = 0.35,
        eta_min: float = 0.80,
        eta_max: float = 0.94,
) -> list[float]:
    """
    Распределяет общий КПД узла по ступеням по форме нормального закона.

    eta_total   - средний КПД узла
    stage_count - число ступеней
    spread      - амплитуда изменения КПД между средними и крайними ступенями
    sigma_rel   - относительная ширина нормального распределения
    eta_min     - нижнее ограничение КПД ступени
    eta_max     - верхнее ограничение КПД ступени

    Возвращает список КПД ступеней длиной stage_count.
    Среднее значение списка равно eta_total с учётом ограничений.
    """

    if stage_count <= 0:
        raise ValueError("stage_count должен быть больше нуля")

    if not (0 < eta_total < 1):
        raise ValueError("eta_total должен быть в диапазоне (0, 1)")

    if stage_count == 1:
        return [eta_total]

    center = (stage_count - 1) / 2
    sigma = sigma_rel * stage_count

    weights = []
    for i in range(stage_count):
        x = (i - center) / sigma
        weights.append(math.exp(-0.5 * x ** 2))

    mean_weight = sum(weights) / stage_count

    eta_stage = [
        eta_total + spread * (w - mean_weight)
        for w in weights
    ]

    eta_stage = [
        min(max(eta, eta_min), eta_max)
        for eta in eta_stage
    ]

    mean_eta = sum(eta_stage) / stage_count
    correction = eta_total - mean_eta

    eta_stage = [
        min(max(eta + correction, eta_min), eta_max)
        for eta in eta_stage
    ]

    return eta_stage


import math


def distribute_stage_work_coefficients(
        L_total: float,
        stage_count: int,
        *,
        u: float | list[float],
        base_coefficients: list[float] | None = None,
        sigma_rel: float = 0.35,
) -> list[float]:
    """
    Распределяет коэффициенты работы ступеней L_rel так, чтобы:

        L_stage[i] = L_rel[i] * u[i]**2
        sum(L_stage) = L_total

    L_total           - суммарная работа узла, Дж/кг
    stage_count       - количество ступеней
    u                 - окружная скорость, м/с:
                        либо одно число, либо список по ступеням
    base_coefficients - исходный профиль L_rel.
                        Если не задан, используется нормальный профиль.
    sigma_rel         - ширина нормального распределения,
                        если base_coefficients не задан.
    """

    if stage_count <= 0:
        raise ValueError("stage_count должен быть больше нуля")

    if L_total <= 0:
        raise ValueError("L_total должна быть больше нуля")

    if isinstance(u, (int, float)):
        u_stage = [float(u)] * stage_count
    else:
        if len(u) != stage_count:
            raise ValueError(
                f"Длина массива u должна быть равна {stage_count}, "
                f"получено {len(u)}."
            )
        u_stage = [float(value) for value in u]

    if base_coefficients is None:
        center = (stage_count - 1) / 2
        sigma = sigma_rel * stage_count

        base_coefficients = []
        for i in range(stage_count):
            x = (i - center) / sigma
            base_coefficients.append(math.exp(-0.5 * x ** 2))
    else:
        if len(base_coefficients) != stage_count:
            raise ValueError(
                f"Длина массива base_coefficients должна быть равна {stage_count}, "
                f"получено {len(base_coefficients)}."
            )

    work_base = sum(
        L_rel * u_i ** 2
        for L_rel, u_i in zip(base_coefficients, u_stage)
    )

    if work_base <= 0:
        raise ValueError("Базовая сумма работ должна быть больше нуля")

    scale = L_total / work_base

    return [
        L_rel * scale
        for L_rel in base_coefficients
    ]