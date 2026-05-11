import math


def velocity_critical(k, R, T):
    return math.sqrt(2 * k / (k + 1) * R * T)

def local_speed_velocity(a_crit, k, lam):
    tau_lam = tau_lambda(lam=lam, k=k)
    return a_crit * math.sqrt((k + 1) / 2 * tau_lam)

def polotrop(pi, T_in, T_out):
    coef = math.log10(pi) / math.log10(T_out / T_in)
    return coef / (coef - 1)

def tau_lambda(lam, k):
    return 1 - (k - 1) / (k + 1) * lam ** 2

def pi_lambda(lam, k):
    pi_lambda = (1 - (k - 1) / (k + 1) * lam**2) ** (k / (k-1))
    return pi_lambda

def lambda_pi( pi_l, k):
    lam = math.sqrt(((k+1)/(k-1)) * (1 - pi_l ** ((k-1)/k)))
    return lam

def z_lambda( lam):
    return 0.5 * (lam + (1 / lam))

def q_lambda(lam: float, k: float) -> float:
    """
    q(λ) = λ * (1 - (k-1)/(k+1) * λ^2)^(1/(k-1)) * ((k+1)/2)^(1/(k-1))
    """
    a = (k - 1.0) / (k + 1.0)
    base = 1.0 - a * lam * lam
    if base <= 0.0:
        return float("nan")
    return lam * (base ** (1.0 / (k - 1.0))) * (((k + 1.0) / 2.0) ** (1.0 / (k - 1.0)))

def lambda_from_q(q: float, k: float, branch: str = "sub", *, tol: float = 1e-6, max_iter: int = 100) -> float:
    """
    Численное обращение q(λ) -> λ методом бисекции (надёжно, без производных).
    branch:
      - "sub" : дозвуковая ветвь, λ ∈ (0, 1]
      - "sup" : сверхзвуковая ветвь, λ ∈ [1, λ_max)
    Возвращается λ > 0.
    """
    if k <= 1.0:
        raise ValueError("k must be > 1")
    if q <= 0.0:
        raise ValueError("q must be > 0")
    # верхняя граница по условию (1 - a*λ^2) > 0
    lam_max = math.sqrt((k + 1.0) / (k - 1.0))
    eps = 1e-12
    q_crit = q_lambda(1.0, k)  # обычно максимум (критика)
    if not math.isfinite(q_crit):
        raise RuntimeError("q(1) is not finite; check formula / k")
    # физический диапазон: q обычно в (0, q_crit]
    if q > q_crit * (1.0 + 1e-12):
        raise ValueError(f"q={q} exceeds q_crit={q_crit} for this definition; no real λ")
        # if branch.lower() in ("sub", "subsonic", "low"):
        #     return 1.0 - 1e-8
        # else:
        #     return 1.0 + 1e-8
    if branch.lower() in ("sub", "subsonic", "low"):
        lo, hi = eps, 1.0
    elif branch.lower() in ("sup", "supersonic", "high"):
        lo, hi = 1.0, lam_max * (1.0 - 1e-12)
    else:
        raise ValueError("branch must be 'sub' or 'sup'")
    # функция для корня
    def f(lam: float) -> float:
        return q_lambda(lam, k) - q
    f_lo, f_hi = f(lo), f(hi)
    # На каждой ветви q(λ) монотонна, поэтому корень единственный.
    # Но на границах возможны численные эффекты — подстраховка.
    if not (math.isfinite(f_lo) and math.isfinite(f_hi)):
        raise RuntimeError("Non-finite values on bracket; adjust eps or check inputs")
    # для "sub": f(lo)<0, f(hi)>=0 ; для "sup": f(lo)>=0, f(hi)<0
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    if f_lo * f_hi > 0.0:
        raise RuntimeError(
            f"Root is not bracketed on {branch} branch: f(lo)={f_lo}, f(hi)={f_hi}. "
            "Check q range or branch selection."
        )
    # бисекция
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if abs(f_mid) <= tol:
            return mid
        # сужение интервала
        if f_lo * f_mid <= 0.0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
        # критерий по длине интервала
        if (hi - lo) <= tol * max(1.0, abs(mid)):
            return 0.5 * (lo + hi)
    raise RuntimeError("Max iterations exceeded")
