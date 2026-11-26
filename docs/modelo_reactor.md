# Ampliación 2: simulador de reactor (calor, agua, O₂ y cinética)

> **Borrador para revisión de Vicente. Aún no hay código.** Este documento extrae el modelo de la literatura antes de implementarlo, para que lo apruebes. Las partes que faltan están marcadas como **[FALTA]**.

## Por qué un modelo publicado y no tus 159 series

El Compost Dataset no informa caudal de aire, O₂, volumen del reactor ni temperatura ambiente. Un balance de energía, de agua y de O₂ no se puede calibrar ni validar con esos datos. Por eso se reproduce un modelo publicado que trae sus propios experimentos de validación. Tu cinética entra como punto de comparación: la corrección de temperatura del modelo es justamente Arrhenius con temperatura de referencia, como en la versión corregida del repo.

## Fuente principal

Papraćanin, E. (2022). *Kinetic parameters estimation and model evaluation for the aerobic bioconversion process of the organic fraction of municipal solid waste.* Kemija u industriji 71(3-4), 185–199. doi:10.15255/KUI.2021.048. Licencia **CC BY 4.0**.

El artículo remite a dos referencias para los balances completos: Papraćanin y Petric (2016), *Bull. Chem. Technol. Bosnia Herzegovina* 47, 39–48 (acceso abierto), y Petric y Selimbašić (2008), *Chem. Eng. J.* 139(2), 304–317, doi:10.1016/j.cej.2007.08.017 (de pago).

## Sistema

Es un reactor piloto de 57 L, aislado, con aire forzado de 0,4 L·min⁻¹·kg_OM⁻¹ y aire ambiente a 22,5 ± 2,5 °C, durante 360 h. Supuestos del autor: volumen de gas y presión constantes, gas saturado de vapor de agua, caudal de aire constante, fases sólido-líquida y gaseosa a la misma temperatura, y sustrato homogéneo al inicio.

## Ecuaciones disponibles en la fuente principal (n.º de ecuación del artículo)

- (1) Degradación: dm_OM/dt = −k · m_OM^n, con m_OM en kg y t en h.
- (2) k = k_T · k_O₂ · k_pH · k_H₂O · k_FAS.
- (3) Temperatura: k_T = α · exp[β (1/293 − 1/T)], con β = Eₐ/R y T en K.
- (4) Oxígeno: k_O₂ = O₂ / [k_O₂(0) · (K_O₂ + O₂)], con O₂ en % v/v.
- (6) Humedad: k_H₂O = 1 / [exp(−γ (1 − S_m) + Δ) + 1], donde S_m es la fracción de materia seca.
- (8) Espacio libre de aire: k_FAS = 1 / [exp(−χ · FAS + ω) + 1].
- (9) Estequiometría: C_aH_bO_cN_d + (4a + b − 2c − 3d)/4 O₂ → a CO₂ + (b − 3d)/2 H₂O + d NH₃.
- (10) CO₂ disuelto: dm_CO₂/dt = Y_CO₂ · dm_OM/dt − k_L a_CO₂ (He_CO₂ f_CO₂ X_CO₂ − p_CO₂).
- (11) Energía de la fase sólido-líquida: dT/dt = [Q_cw + Q_G − Q_c^T − Σ entalpías de los gases transferidos − (T − 273,15)(c_pw dm_w/dt + c_pOM dm_OM/dt)] / (c_pw m_w + c_pOM m_OM + c_pNT m_NT).
- Tabla 3: balances de gases disueltos (CO₂, O₂, NH₃), de agua en el material, de cinco gases en la fase gaseosa y de temperatura del gas. En total son 12 variables de estado.

## Parámetros ajustados por el autor (Tabla 5, mezcla M1)

| Parámetro | Valor | Unidad | IC 95 % |
| --- | ---: | --- | --- |
| α | 1,2688·10⁻⁵ | según n | 1,12–1,42·10⁻⁵ |
| β = Eₐ/R | 4544,1 | K (Eₐ = 37,8 kJ/mol) | 4401–4679 |
| n | ≈ 1,9 | — | 1,65–2,10 |
| Y_CO₂ | 3,17 | kg CO₂/kg | 2,96–3,36 |
| Q_c | 1753,3 | J/kg_OM | 1704–1796 |
| k_O₂(0) | 0,1554 | — | 0,129–0,182 |
| K_O₂ | 0,868 | % v/v | 0,37–1,36 |
| γ, Δ | 17,80; 6,86 | — | — |
| χ, ω | 23,67; 3,50 | — | — |

Condiciones iniciales de M1: m_OM = 6,3 kg, agua = 11,35 kg, O₂ = 5,84·10⁻⁵ kmol, CO₂ = 5·10⁻⁶ kmol y T = 295 K.

## [FALTA] para implementar

Estos datos están en Papraćanin y Petric (2016) o en Petric y Selimbašić (2008), no en la fuente principal:

1. Q_G: el calor de reacción por kg de materia orgánica degradada.
2. Q_cw: la forma de la pérdida por paredes, con su coeficiente y su área.
3. c_pw, c_pOM, c_pNT y los c_p de los gases.
4. k_L a, He y f para O₂, CO₂ y NH₃, la constante Y_w y la correlación de presión de vapor de la Tabla 3.
5. k_pH, que el artículo nombra pero no define.
6. FAS inicial, composición elemental (a, b, c, d) y masa inorgánica m_NT de cada mezcla.
7. Una incoherencia de unidades: α figura en kg⁻⁰·⁵ h⁰·⁵, pero con n ≈ 1,9 la ecuación (1) exige kg⁻⁰·⁹ h⁻¹.

## Caso de validación propuesto

Simular las mezclas M2 y M3 con los parámetros de M1 y comparar con las Figuras 2 y 3 del artículo, digitalizadas. Las referencias son las del autor: RMSE de temperatura de 2,89 K (M2) y 2,21 K (M3), y RMSE de materia orgánica de 0,20 y 0,21 kg. El simulador queda validado si llega a esos órdenes de magnitud.

## Qué aportaría al diseño

Con el modelo validado se puede estudiar la relación entre caudal de aire y temperatura máxima (higienización sobre 55 °C), la pérdida de humedad, que puede frenar la degradación, y el tiempo para llegar a una conversión dada. Eso responde la pregunta de la memoria, el tiempo de residencia, con física y no solo con un ajuste empírico.
