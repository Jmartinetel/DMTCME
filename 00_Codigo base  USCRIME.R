library(tidyverse)
library(janitor)
library(car)

library(glmnet)
library(MASS)


# Cargar datos
data(UScrime)
df <- UScrime %>%
  dplyr::select(-Po1)


set.seed(197)
df <- df %>% slice_sample(n = 45)
# Variable respuesta
y <- df$y

# Matriz de regresores
X <- df %>% dplyr::select(-y)


# Correlaciones entre X
cor(X)

# Modelo OLS
modelo <- lm(y ~ ., data = df)


summary(modelo)
# VIF para multicolinealidad
car::vif(modelo)


X <- df %>% dplyr::select(-y) %>% as.matrix()

Xc <- scale(X, center = TRUE, scale = FALSE)
Yc <- y - mean(y)

n <- nrow(Xc)
k <- ncol(Xc)




# svd ---------------------------------------------------------------------

svd_X <- svd(Xc, nu = k, nv = k)

U <- svd_X$u
d <- svd_X$d

V <- svd_X$v




# sigma -------------------------------------------------------------------



P_res <- diag(n) - U %*% t(U)

sigma2_hat <- drop(t(Yc) %*% P_res %*% Yc) / (n - k - 1)




# toro y wllace -----------------------------------------------------------

prueba <- readxl::read_xlsx("data/tabla_cme_final.xlsx") %>%
  clean_names()

resultados_test <- tibble(
  l = integer(),
  s = integer(),
  T_stat = numeric(),
  T_alpha = numeric(),
  decision = character()
)

ultimo_l_aceptado <- 0


for (l in 1:(k-1)) {
  # browser()
  
  # Componentes a eliminar = últimos l
  U_k_s <- U[, (k - l + 1):k, drop = FALSE]
  
  T_stat <- drop(
    t(Yc) %*% U_k_s %*% t(U_k_s) %*% Yc
  ) / (l * sigma2_hat)
  
  T_alpha <- prueba %>%
    filter(n_k == (n - k - 1)) %>%
    pull(!!sym(paste0("x", l)))
  
  decision <- ifelse(T_stat < T_alpha, "Aceptar H0", "Rechazar H0")
  
  if (decision == "Aceptar H0") {
    ultimo_l_aceptado <- l
  }
  
  resultados_test <- add_row(
    resultados_test,
    l = l,
    s = k - l,
    T_stat = T_stat,
    T_alpha = T_alpha,
    decision = decision
  )
}

resultados_test


s_sel <- k - ultimo_l_aceptado
s_sel




# eliinar -----------------------------------------------------------------

# Subespacio reducido
V_s <- V[, 1:s_sel, drop = FALSE]
D_s <- diag(svd_X$d[1:s_sel])
D_s_inv <- solve(D_s)

# Regresión en el espacio reducido
X_reducido <- Xc %*% V_s
modelo_dtmce <- lm(y ~ X_reducido)

gamma_hat <- coef(modelo_dtmce)[-1]


# Coeficientes en la base original
beta_dtmce_completo <- V_s %*% gamma_hat



calcular_ur <- function(V_s, gamma_hat, D_s_inv, r, variable_names) {
  
  combs <- combn(seq_len(nrow(V_s)), r)
  
  ur_df <- tibble(
    pair = character(ncol(combs)),
    ur_value = numeric(ncol(combs))
  )
  
  for (i in seq_len(ncol(combs))) {
    
    idx <- combs[, i]
    V_r <- V_s[idx, , drop = FALSE]
    
    ur <- drop(
      t(gamma_hat) %*% t(V_r) %*%
        (V_r %*% D_s_inv %*% t(V_r)) %*%
        V_r %*% gamma_hat
    )
    
    ur_df$pair[i] <- paste(variable_names[idx], collapse = "-")
    ur_df$ur_value[i] <- ur
  }
  
  arrange(ur_df, ur_value)
}

# Regresión en el espacio reducido
X_reducido <- Xc %*% V_s

modelo_dtmce <- lm(y ~ X_reducido)

gamma_hat <- coef(modelo_dtmce)[-1]

# Coeficientes en la base original (estimador restringido)
beta_dtmce <- V_s %*% gamma_hat


sigma2_dtmce <- sum(residuals(modelo_dtmce)^2) /
  (n - s_sel - 1)

tr_ECM_dtmce <- sigma2_dtmce * sum(1 / svd_X$d[1:s_sel]^2)

# mansflied ---------------------------------------------------------------


D_s_inv <- solve(D_s)


variable_names <- colnames(Xc)
temp_V_s <- V_s

mansfield_results <- list()

max_r <- nrow(temp_V_s) - 1

for (r in 1:max_r) {
  # browser()
  
  if (nrow(temp_V_s) < r) break
  
  ur_table <- calcular_ur(
    V_s = temp_V_s,
    gamma_hat = gamma_hat,
    D_s_inv = D_s_inv,
    r = 1,
    variable_names = variable_names
  )
  
  mansfield_results[[r]] <- ur_table
  
  # Variable con menor contribución
  var_to_remove <- ur_table$pair[1]
  idx_remove <- match(var_to_remove, variable_names)
  
  # Eliminar
  temp_V_s <- temp_V_s[-idx_remove, , drop = FALSE]
  variable_names <- variable_names[-idx_remove]
  
  cat("Iteración", r, "- Eliminando:", var_to_remove, "\n")
  cat("Variables restantes:", variable_names, "\n\n")
}

X_final <- Xc[, variable_names, drop = FALSE]

modelo_final_2 <- lm(y ~ ., data = as.data.frame(X_final))
summary(modelo_final)
summary(modelo_final_2)

car::vif(modelo_final)



# ============================================================
# COMPARACIÓN DE MODELOS
# ============================================================

library(car)
library(glmnet)
library(tibble)
library(dplyr)

# ============================================================
# OLS COMPLETO
# ============================================================

df_ols <- as.data.frame(Xc)
df_ols$y <- y

modelo_ols <- lm(y ~ ., data = df_ols)

# VIF
VIF_ols <- car::vif(modelo_ols)
VIF_ols_max <- max(VIF_ols)

# CME empírico
y_hat_ols <- predict(modelo_ols)
CME_ols <- mean((y - y_hat_ols)^2)

# Traza del ECM (OLS)
svd_X <- svd(Xc)
sigma2_ols <- sum(residuals(modelo_ols)^2) / df.residual(modelo_ols)
tr_ECM_ols <- sigma2_ols * sum(1 / svd_X$d^2)

R2_ols <- summary(modelo_ols)$adj.r.squared


# ============================================================
# DTMCE (GERIG–ZÁRATE + MANSFIELD)
# ============================================================

# CME empírico (modelo final en variables)
y_hat_dtmce <- predict(modelo_final)
CME_dtmce <- mean((y - y_hat_dtmce)^2)

# VIF post-Mansfield
VIF_dtmce <- car::vif(modelo_final)
VIF_dtmce_max <- max(VIF_dtmce)

# Traza del ECM (MODELO EN COMPONENTES, no en variables)
sigma2_dtmce <- sum(residuals(modelo_dtmce)^2) / (n - s_sel - 1)
tr_ECM_dtmce <- sigma2_dtmce * sum(1 / svd_X$d[1:s_sel]^2)

R2_dtmce <- summary(modelo_final)$adj.r.squared


# ============================================================
# RIDGE
# ============================================================

set.seed(123)

cv_ridge <- cv.glmnet(
  x = as.matrix(Xc),
  y = y,
  alpha = 0,
  standardize = FALSE
)

beta_ridge <- as.vector(coef(cv_ridge, s = "lambda.min"))[-1]
names(beta_ridge) <- colnames(Xc)

y_hat_ridge <- as.vector(predict(cv_ridge, as.matrix(Xc), s = "lambda.min"))
CME_ridge <- mean((y - y_hat_ridge)^2)

R2_ridge <- 1 - sum((y - y_hat_ridge)^2) / sum((y - mean(y))^2)


# ============================================================
# LASSO
# ============================================================

cv_lasso <- cv.glmnet(
  x = as.matrix(Xc),
  y = y,
  alpha = 1,
  standardize = FALSE
)

beta_lasso <- as.vector(coef(cv_lasso, s = "lambda.min"))[-1]
names(beta_lasso) <- colnames(Xc)

y_hat_lasso <- as.vector(predict(cv_lasso, as.matrix(Xc), s = "lambda.min"))
CME_lasso <- mean((y - y_hat_lasso)^2)

n_vars_lasso <- sum(beta_lasso != 0)

R2_lasso <- 1 - sum((y - y_hat_lasso)^2) / sum((y - mean(y))^2)




# ============================================================
# TABLA COMPARATIVA DE MÉTRICAS
# ============================================================

tabla_comparativa <- tibble(
  Metodo = c("OLS", "DTMCE", "Ridge", "Lasso"),
  R2_Ajustado = c(R2_ols, R2_dtmce, R2_ridge, R2_lasso),
  CME = c(CME_ols, CME_dtmce, CME_ridge, CME_lasso),
  VIF_max = c(VIF_ols_max, VIF_dtmce_max, NA, NA),
  Traza_ECM = c(tr_ECM_ols, tr_ECM_dtmce, NA, NA),
  Comentario = c(
    "OLS completo",
    paste0("s = ", s_sel, " componentes"),
    "Penalización L2",
    paste0("L1 (", n_vars_lasso, " vars)")
  )
)

tabla_comparativa

car::vif(modelo)
# ============================================================
# TABLA DE COEFICIENTES
# ============================================================

beta_ols <- coef(modelo_ols)[-1]
beta_dtmce_final <- coef(modelo_final)[-1]
beta_masnifled_ols <- coef(modelo_final_2)[-1]
vars <- colnames(Xc)

tabla_coeficientes <- tibble(
  Variable = vars,
  OLS = beta_ols[vars],
  DTMCE = beta_dtmce[,1],
  Masnfield =  ifelse(vars %in% names(beta_masnifled_ols),
                      beta_masnifled_ols[vars], 0),
  DTMCE_MANSFILED = ifelse(vars %in% names(beta_dtmce_final),
                           beta_dtmce_final[vars], 0),
  Ridge = beta_ridge[vars],
  Lasso = beta_lasso[vars]
)

tabla_coeficientes

