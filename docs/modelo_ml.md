# Ampliación 1: predecir k sin ensayo cinético (aprendizaje automático)

> **Borrador para revisión de Vicente.** Este documento amplía lo que la memoria dejó fuera de alcance: modelos de aprendizaje automático y predicción fuera del rango de los datos. Las elecciones de método están marcadas con **[Propuesta]** y esperan tu visto bueno.

## Pregunta de ingeniería

Al diseñar una planta, ¿basta con conocer el sustrato y las condiciones iniciales para estimar k? ¿O hace falta un ensayo cinético? Y si alguien pide una predicción para condiciones que el modelo nunca vio, ¿cómo se le avisa?

## Datos

Hay 150 de las 159 series. Se descartan las que tienen un ajuste de primer orden con R² < 0,5 **[Propuesta]** o k ≤ 10⁻⁴ d⁻¹, porque su k no representa a la serie. El k es el de la versión corregida (`results/ajustes_por_serie.csv`).

| Variable | Tipo | Unidad | Rango en los datos |
| --- | --- | --- | --- |
| Sustrato (cocina, estiércol, lodos) | entrada, categórica | — | 45 / 47 / 58 series |
| pH inicial | entrada | — | 5,1–8,8 |
| C/N inicial | entrada | — | 8,5–54 |
| Humedad inicial | entrada | % | 48–75 |
| Temperatura media del ensayo | entrada | °C | 18–61 |
| log₁₀ k | salida | k en d⁻¹ | k entre 0,003 y 0,6 d⁻¹ |

**[Propuesta]** Se usa la temperatura media aunque sea una variable del proceso y no una condición inicial. La razón es que Arrhenius muestra que la temperatura importa. Si se quita, el modelo sirve antes de operar, pero predice peor.

## Modelos

Hay tres modelos **[Propuesta]**: la mediana de entrenamiento como referencia, una regresión Ridge y un bosque aleatorio (300 árboles, hoja mínima de 5). Ninguno se ajusta con más detalle, porque con 150 datos eso sería sobreajustar la validación.

## Validación

Se usan dos esquemas, porque responden preguntas distintas:

1. **Series nuevas.** Es una validación cruzada de 5 pliegues por grupo. Las series con datos (t, C) idénticos (los duplicados de `supuestos.md`, L5) caen siempre en el mismo pliegue, para que el modelo no se evalúe con datos que ya vio.
2. **Sustrato nuevo.** Se entrena con dos sustratos y se predice el tercero. Esto mide la predicción fuera de rango.

Métricas: RMSE en log₁₀ k, R² y fracción de predicciones dentro de un factor 2 del k ajustado.

## Resultado

| Validación | Modelo | RMSE (log₁₀) | R² | Dentro de factor 2 |
| --- | --- | ---: | ---: | ---: |
| Series nuevas | Mediana | 0,50 | −0,03 | 45 % |
| Series nuevas | Ridge | 0,50 | −0,03 | 49 % |
| Series nuevas | Bosque aleatorio | 0,45 | 0,17 | 60 % |
| Sustrato nuevo | Mediana | 0,49 | −0,02 | 45 % |
| Sustrato nuevo | Ridge | 0,52 | −0,11 | 44 % |
| Sustrato nuevo | Bosque aleatorio | 0,54 | −0,21 | 54 % |

**Lectura.** Las condiciones iniciales explican poco de k. El bosque aleatorio explica un 17 % de la varianza en series nuevas y acierta dentro de un factor 2 en 6 de cada 10. Con un sustrato que no vio, predice peor que la mediana. Para dimensionar un reactor, hay que medir la cinética del sustrato real. Un modelo entrenado con literatura sirve como primera estimación, no como dato de diseño.

## Dominio de aplicabilidad (predicción fuera de rango)

Una entrada queda **dentro** del dominio cuando se cumplen dos condiciones. La primera es que cada variable numérica esté dentro del rango de entrenamiento. La segunda es que su distancia media a los 5 vecinos más cercanos, con las variables estandarizadas, no supere el percentil 95 de esa misma distancia en los datos de entrenamiento **[Propuesta]**. Fuera del dominio, la app muestra la predicción con una advertencia.

En validación cruzada, las 19 series que cayeron fuera del dominio tuvieron un error mediano de 0,42 en log₁₀, es decir, un factor 2,6. Las 131 que quedaron dentro tuvieron 0,20, un factor 1,6. Dentro de un factor 2 quedaron el 37 % de las de afuera y el 63 % de las de adentro. El aviso separa las predicciones menos confiables.

## Supuestos y límites

- k de primer orden es un buen resumen de la serie, por eso se filtra por R².
- Las series de la literatura representan las condiciones en que se usará el modelo.
- No hay variables de aireación, tamaño de partícula ni tipo de reactor en el dataset. Probablemente explican parte de la varianza que el modelo no captura.

Código: `src/compostaje/ml.py`. Tests: `tests/test_ml.py`. Resultados: `results/ml_*.csv` y `docs/figuras/ml_paridad.png`.
