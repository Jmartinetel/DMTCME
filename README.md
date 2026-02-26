# Selección de variables mediante el método DTMCE

Este repositorio contiene la implementación en R del algoritmo utilizado en la tesis para la selección de variables en modelos de regresión lineal con problemas de multicolinealidad, basado en el enfoque de estimadores direccionalmente minimax en la traza del error cuadrático medio (DTMCE).

El procedimiento combina la reducción de componentes principales con el criterio de selección de Mansfield, siguiendo la metodología desarrollada por Gerig y Zárate.

## Implementación principal

La función principal se encuentra en:

`lib/gerig_zarate_algoritmo.R`

y corresponde a:

`fit_dtmce_mansfield()`



Esta función ajusta un modelo lineal con selección de variables basada en el método DTMCE, identificando y eliminando aquellas variables asociadas a componentes principales redundantes o inestables.

## Archivos incluidos

- **00_Codigo base USCRIME.R**  
  Código utilizado en el ejemplo principal presentado en la tesis.

- **01_ejemplo_uscrime.R**  
  Aplicación del método a la base de datos USCRIME.

- **02_ejemplo_Boston.R**  
  Aplicación del método a la base Boston Housing.

- **03_ejemplos_mtcars.R**  
  Aplicación del método a la base mtcars.

## Objetivo

Este repositorio tiene como propósito:

- Proporcionar una implementación reproducible del algoritmo desarrollado
- Facilitar su aplicación a otros conjuntos de datos
- Servir como complemento computacional del trabajo de tesis
