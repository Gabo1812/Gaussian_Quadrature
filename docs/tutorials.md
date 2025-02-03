# Tutorials

A continuación se muestra cómo utilizar la implementación de la cuadratura de Gauss-Legendre para calcular una integral definida utilizando Python, para distintos valores de N.

### Paso 1: Definir la función a integrar

Definimos la función que queremos integrar. En este caso, utilizaremos la función:

\[
f(x) = x^6 - \sin(2x) \cdot x^2
\]

```python
def func_integrar(x):
    return x**6 - np.sin(2*x) * x**2
```
### Paso 2: Preparar el Código Principal

En el siguiente paso, definimos el intervalo de integración en el que vamos a calcular la integral. En nuestro caso, el intervalo será de [1,3].
```python
lim_inf = 1  # Límite inferior de la integral
lim_sup = 3  # Límite superior de la integral
```
Luego, utilizamos la función ```gaussxw(N)```  para calcular los puntos de colocación y los pesos asociados a la cuadratura de Gauss, donde N es el número de puntos de colocación que utilizamos. Posteriormente, reescalamos estos puntos y pesos al intervalo deseado [1,3] utilizando la función ```gaussxwab(a, b, x, w)```.

Finalmente, evaluamos la integral para distintos valores de N, sumando los productos de los valores de la función en los puntos de integración escalados por sus respectivos pesos.

```python
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
```

Con estos pasos, obtendremos un valor aproximado para la integral de la función 
f(x) en el intervalo [1, 3] usando la cuadratura de Gauss-Legendre con diferentes valores de N.


En este ejemplo el output sería:
```
Para N=2:
Puntos de muestreo: [0.57735027 0.57735027]
Pesos: [1. 1.] 

Para N=2, escalado:
Puntos escalados: [2.57735027 2.57735027]
Pesos escalados: [1. 1.] 

Valor de la integral para N=2: 598.24

Para N=3:
Puntos de muestreo: [ 7.74596669e-01 -1.57179011e-16 -3.39183690e-17]
Pesos: [0.55555556 0.88888889 0.88888889] 

Para N=3, escalado:
Puntos escalados: [2.77459667 2.         2.        ]
Pesos escalados: [0.55555556 0.88888889 0.88888889] 

Valor de la integral para N=3: 375.494

Para N=4:
Puntos de muestreo: [ 0.86113631  0.33998104 -0.86113631 -0.33998104]
Pesos: [0.34785485 0.65214515 0.34785485 0.65214515] 

Para N=4, escalado:
Puntos escalados: [2.86113631 2.33998104 1.13886369 1.66001896]
Pesos escalados: [0.34785485 0.65214515 0.34785485 0.65214515] 

Valor de la integral para N=4: 317.345

Para N=5:
Puntos de muestreo: [ 9.06179846e-01  5.38469310e-01 -8.24468998e-17  5.38469310e-01
 -5.38469310e-01]
Pesos: [0.23692689 0.47862867 0.56888889 0.47862867 0.47862867] 

Para N=5, escalado:
Puntos escalados: [2.90617985 2.53846931 2.         2.53846931 1.46153069]
Pesos escalados: [0.23692689 0.47862867 0.56888889 0.47862867 0.47862867] 

Valor de la integral para N=5: 448.116

Para N=6:
Puntos de muestreo: [ 0.93246951  0.66120939  0.23861919 -0.66120939  0.93246951 -0.66120939]
Pesos: [0.17132449 0.36076157 0.46791393 0.36076157 0.17132449 0.36076157] 

Para N=6, escalado:
Puntos escalados: [2.93246951 2.66120939 2.23861919 1.33879061 2.93246951 1.33879061]
Pesos escalados: [0.17132449 0.36076157 0.46791393 0.36076157 0.17132449 0.36076157] 

Valor de la integral para N=6: 414.077

Para N=7:
Puntos de muestreo: [ 9.49107912e-01  7.41531186e-01  4.05845151e-01  7.54442507e-33
 -9.49107912e-01 -4.05845151e-01 -7.41531186e-01]
Pesos: [0.12948497 0.27970539 0.38183005 0.41795918 0.12948497 0.38183005
 0.27970539] 

Para N=7, escalado:
Puntos escalados: [2.94910791 2.74153119 2.40584515 2.         1.05089209 1.59415485
 1.25846881]
Pesos escalados: [0.12948497 0.27970539 0.38183005 0.41795918 0.12948497 0.38183005
 0.27970539] 

Valor de la integral para N=7: 317.344
```

Podemos observar que para $𝑁=4$ obtenemos el valor correcto de la integral, ya que, según la teoría de la cuadratura de Gauss-Legendre, esta es exacta para polinomios de grado $(2𝑁−1)$.

Es decir, con $𝑁=4$, la cuadratura es exacta para polinomios de grado 7 o inferior.