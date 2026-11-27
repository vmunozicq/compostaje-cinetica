# Ampliación 3: dimensionamiento preliminar de una planta

> **Especificación aprobada por Vicente el 2026-11-26** (S-D1 a S-D5). Esta ampliación no altera ninguna cifra existente del repositorio.

## Pregunta de ingeniería

Tengo X toneladas por día de un residuo. ¿Cuántos días tiene que estar en compostaje, cuánto volumen y área ocupa, cuánto aire hay que soplar y cuánta agua se pierde? Es la pregunta de la memoria, el tiempo de residencia, convertida en un número de planta con su incertidumbre.

## Ecuaciones

1. **Tiempo de residencia** (primer orden, `docs/modelo.md`). De C(t) = C∞ + (C₀ − C∞)·e^(−kt), la fracción X de lo degradable se alcanza en t = −ln(1 − X)/k.
2. **k del sustrato.** Se usan los percentiles 25, 50 y 75 de k de primer orden en las series del sustrato con R² ≥ 0,5. Es el mismo filtro del modelo ML: 45 series de cocina, 47 de estiércol y 58 de lodos. El P25 es el escenario conservador de diseño. El usuario puede ingresar en su lugar un k medido para su residuo, por ejemplo en la pestaña "Ajustar mis datos".
3. **Materia orgánica degradada.** MO_deg = Q · (1 − H) · f_MO · f_deg · X, donde f_deg es la mediana de (C₀ − C∞)/C₀ en las mismas series.
4. **Estequiometría** (Petric y Selimbašić, 2008, la misma del simulador): O₂ = 0,705 · MO_deg, CO₂ = 1,276 · MO_deg, calor = 15,4 MJ/kg · MO_deg.
5. **Aire estequiométrico** = O₂ / (0,2315 · 1,204 kg/m³), con aire a 20 °C y 1 atm.
6. **Aire de diseño y evaporación.** Se calculan por kg de MO degradada en el simulador validado (experimento 2, 0,18 m³/h): 47,9 m³ de aire y 2,75 kg de agua evaporada por kg de MO degradada. Se escalan linealmente con MO_deg.
7. **Volumen** = Q · t / ρ_aparente. **Área de pila** = volumen / altura.

## Supuestos

- **S-D1.** f_deg y k se definen sobre la concentración de TOC (% b.s.), no sobre la masa. Es la misma variable que ajusta la memoria. La MO degradada que resulta es una estimación de orden de magnitud.
- **S-D2.** El aire y la evaporación por kg de MO degradada del reactor de 32 L se trasladan a escala de planta. En un reactor de laboratorio las pérdidas por pared pesan más que en una pila grande, así que el aire real de control térmico puede ser distinto. La intensidad de aireación se conserva por kg degradado, no por m³.
- **S-D3.** El volumen usa la masa húmeda de entrada sin descontar la pérdida de masa durante el proceso. Queda del lado conservador.
- **S-D4.** La densidad aparente y la altura de pila no están en las fuentes del repositorio. Son entradas del usuario; 600 kg/m³ y 2 m son solo valores de ejemplo (aprobado 2026-11-26).
- **S-D5.** Los k vienen de ensayos de laboratorio de la literatura, con su temperatura. No se corrigen por la temperatura de la planta.

## Caso de ejemplo (10 t/d de residuo de cocina, humedad 60 %, MO 80 % b.s., X = 90 %)

| Escenario | k (d⁻¹) | Tiempo (d) | Volumen (m³) | Área a 2 m (m²) |
| --- | ---: | ---: | ---: | ---: |
| Conservador (P25) | 0,033 | 69 | 1.148 | 574 |
| Típico (P50) | 0,077 | 30 | 501 | 250 |
| Optimista (P75) | 0,170 | 14 | 225 | 113 |

Para los tres escenarios, la MO degradada es de 1.069 kg/d y el O₂ de 753 kg/d. El aire estequiométrico es de 113 m³/h y el aire de control térmico de unos 2.130 m³/h. Se evaporan cerca de 2,9 t/d de agua.

La dispersión entre P25 y P75 es de 5 veces en el tiempo de residencia. Por eso la conclusión de la memoria se mantiene: para diseñar en serio hay que medir el k del residuo real.

## Fuera de alcance

Volteos, lixiviados, emisiones de NH₃ y olores, costos, y cualquier límite normativo (DS u otra norma chilena), que no se evalúa aquí.
