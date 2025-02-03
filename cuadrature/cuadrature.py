import numpy as np

def gaussxw(N):
    """
    Calcula los puntos de colocación (x_k) y los pesos (w_k) para la cuadratura de Gauss.

    Este método utiliza una aproximación inicial seguida de un método de Newton para 
    calcular las raíces de los polinomios de Legendre. Luego, se calculan los pesos 
    asociados a estas raíces.

    Args:
        N (int): El número de puntos de colocación (orden de la cuadratura).

    Returns:
        tuple: Un tuple con dos elementos:
            - x (ndarray): Los puntos de colocación de Gauss.
            - w (ndarray): Los pesos asociados a los puntos de colocación.
    
    """
    # Aproximación inicial de los puntos
    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    # Método de Newton para calcular las raíces de los polinomios de Legendre
    epsilon = 1e-15  # Tolerancia para la convergencia
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, dtype = float)
        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = np.max(np.abs(dx))

    # Cálculo de los pesos asociados a los puntos
    w = 2 * (N + 1) * (N + 1) / (N * N * (1 - x * x) * dp * dp)

    return x, w

def gaussxwab(a, b, x, w):
    """
    Reescala los puntos y pesos de la cuadratura de Gauss al intervalo [a, b].

    Esta función toma los puntos de colocación (x) y los pesos (w) calculados para 
    el intervalo estándar [-1, 1] y los transforma para que correspondan al intervalo 
    [a, b].

    Args:
        a (float): Límite inferior del intervalo.
        b (float): Límite superior del intervalo.
        x (ndarray): Los puntos de colocación de Gauss en el intervalo [-1, 1].
        w (ndarray): Los pesos asociados a los puntos de colocación en el intervalo [-1, 1].

    Returns:
        tuple: Un tuple con dos elementos:
            - x (ndarray): Los puntos de colocación reescalados al intervalo [a, b].
            - w (ndarray): Los pesos reescalados para el intervalo [a, b].
    
    """
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def func_integrar(x):
    """
    Función a integrar en la cuadratura de Gauss.

    Esta es la función definida para realizar la integración utilizando la 
    cuadratura de Gauss.

    Args:
        x (ndarray): Los puntos en los cuales se evalúa la función.

    Returns:
        f(x) (ndarray): Los valores de la función evaluados en los puntos x.

    """
    return x**6 - np.sin(2*x) * x**2

# Main Code
lim_inf = 1  # Límite inferior de la integral
lim_sup = 3  # Límite superior de la integral

# Evaluación para distintos valores de N
for N in [2, 3, 4, 5, 6, 7]:
    print(f"Para N={N}:")
    
    # Calcular puntos de colocación y pesos utilizando la cuadratura de Gauss
    p_muest, peso = gaussxw(N)
    print("Puntos de muestreo:", p_muest)
    print("Pesos:", peso, "\n")

    # Escalar los puntos y pesos al intervalo [1, 3]
    pto_esc, peso_esc = gaussxwab(lim_inf, lim_sup, p_muest, peso)
    print(f"Para N={N}, escalado:")
    print("Puntos escalados:", pto_esc)
    print("Pesos escalados:", peso_esc, "\n")

    # Evaluar la integral usando los puntos y pesos escalados
    integral = np.sum(func_integrar(pto_esc) * peso_esc)
    print(f"Valor de la integral para N={N}: {round(integral, 3)}\n")
