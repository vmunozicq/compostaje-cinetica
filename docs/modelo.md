# Modelo (una página)

> Extraído de mi memoria de título: *Modelación cinética multimodelo del compostaje aeróbico de sólidos a partir de literatura: ajuste de parámetros y guía para diseño experimental* (Universidad de Concepción, 2025; profesor guía: Dr. Claudio Zaror Z.). Secciones 3.3–3.6 y Anexos A y B. Donde el texto de la memoria y los scripts que generaron sus resultados no coinciden, este repositorio sigue a los scripts y la diferencia queda registrada en [`supuestos.md`](supuestos.md).

## Pregunta de ingeniería

¿Qué modelo cinético describe mejor la degradación de la materia orgánica de un sustrato dado (residuos de cocina, estiércol o lodos), y cuán inciertos son sus parámetros, para dimensionar el tiempo de residencia de un reactor de compostaje?

## Entradas y salidas

| Variable | Tipo | Unidad | Rango válido |
| --- | --- | --- | --- |
| t | entrada | d | ≥ 0 |
| C(t): TOC u OM | entrada | % base seca | 0–100 |
| T(t) | entrada (solo Arrhenius) | °C | 10–75 |
| k, k₁, k₂ | salida | d⁻¹ | ≥ 0 |
| Rm, λ | salida | %·d⁻¹, d | ≥ 0 |
| A, Eₐ | salida | d⁻¹, J·mol⁻¹ | A ∈ [1, 10⁹], Eₐ ∈ [1, 150] kJ/mol |
| R², RMSE, AIC, AICc, BIC | salida | — | — |

## Ecuaciones

C₀ es el primer dato medido (fijo, no se ajusta). C∞ es la fracción no degradable.

1. **Primer orden** (Saldarriaga et al., 2018; Ibrić et al., 2021): dC/dt = −k (C − C∞) ⇒ C(t) = C∞ + (C₀ − C∞) e^(−k t). Parámetros: k, C∞ (p = 2).
2. **Bifásico** (fracción lábil y recalcitrante en paralelo): C(t) = C∞ + (C₀ − C∞) [f e^(−k₁ t) + (1 − f) e^(−k₂ t)], 0 ≤ f ≤ 1, k₁ ≥ k₂. Parámetros: k₁, k₂, f, C∞ (p = 4).
3. **Gompertz modificado** (Zwietering et al., 1990), aplicado a la degradación acumulada Y = C₀ − C: Y(t) = P exp{−exp[(Rm e / P)(λ − t) + 1]}. Parámetros: P, Rm, λ (p = 3).
4. **Arrhenius** (no isotérmico): k(T) = A exp(−Eₐ / (R T)), dC/dt = −k(T(t)) (C − C∞), C(0) = C₀, R = 8,314 J mol⁻¹ K⁻¹, T(t) interpolada linealmente entre mediciones. Se resuelve con `odeint`. Parámetros: A, Eₐ, C∞ (p = 3).

**Ajuste:** mínimos cuadrados no lineales sin ponderar, SSE = Σ (C_obs − C_pred)². `curve_fit` (Trust Region Reflective con cotas) para 1–3; L-BFGS-B sobre el SSE para 4.

**Métricas (Anexo B):** R² = 1 − SSE/SST; RMSE; RMSE normalizado = RMSE / C̄ × 100 %; AIC = n ln(SSE/n) + 2p; AICc = AIC + 2p(p + 1)/(n − p − 1); BIC = n ln(SSE/n) + p ln n. El mejor modelo es el de menor AICc entre los que tienen R² > 0 y n > p + 1.

**Incertidumbre (sección 3.5):** bootstrap de pares con reemplazo, B = 1000; se reajustan primer orden y Gompertz y se reporta el IC 95 % por percentiles (2,5 y 97,5).

## Supuestos

1. Humedad y oxígeno no limitantes en los experimentos seleccionados, de modo que la cinética depende de t, C y T (sección 1.2).
2. TOC u OM (%) representan el sustrato degradable. Se usa TOC si la serie tiene ≥ 5 datos de TOC; si no, OM.
3. Una serie experimental se identifica por (tipo de residuo, pH, C/N y humedad iniciales); cada vez que el tiempo retrocede comienza una réplica nueva.
4. Filtro de calidad: puntos con T ∈ [10, 75] °C, humedad ∈ [30, 75] %, pH ∈ [5,5; 9] y C/N ∈ [8, 40]; series con ≥ 5 puntos completos (t, C, T) y variación ≥ 2 puntos porcentuales.
5. Errores homocedásticos e independientes (pesos wᵢ = 1).

## Caso de validación

Reproducir con el dataset fijado los resultados de la memoria (tablas originales en `data/reference/`):

- 417 series tras el filtro físico, 3.362 puntos y 159 series ajustables.
- R² de primer orden, bifásico y Gompertz idénticos en las 159 series (tolerancia 10⁻⁴; mismo algoritmo, cambian solo las versiones de SciPy).
- Ganadores por AICc de la memoria (primer orden 71, Gompertz 69, bifásico 18, Arrhenius 1) al usar el mismo cálculo de N y SST que el script original.

Tests: [`tests/test_validacion.py`](../tests/test_validacion.py).

## Fuera de alcance

Transferencia de calor y masa, aireación, balance de agua, cinética microbiana explícita, predicción fuera del rango de los datos y modelos de aprendizaje automático (discutidos como trabajo futuro en la memoria).
