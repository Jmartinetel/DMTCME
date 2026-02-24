# ============================================================
# USCRIMe
# ============================================================

library(MASS)
library(dplyr)
library(readxl)
library(janitor)
library(car)

source("lib/gerig_zarate_algoritmo.R")

data(Boston)

df <- Boston

y <- df$medv
X <- df %>% select(-medv)



# ============================================================
# MODELO COMPLETO
# ============================================================

modelo_full <- lm(y ~ ., data = df)

print(vif(modelo_full))


# ============================================================
# ALGORITMO GZ + MANSFIELD
# ============================================================

resultado <- fit_dtmce_mansfield(
  X = X,
  y = y,
  alpha = 0.05
)

vars_final <- resultado$vars_final

# ============================================================
# MODELO REDUCIDO
# ============================================================

modelo_reducido <-resultado$modelo_mansfield

cat("\n========== VIF MODELO REDUCIDO ==========\n")
print(vif(modelo_reducido))
