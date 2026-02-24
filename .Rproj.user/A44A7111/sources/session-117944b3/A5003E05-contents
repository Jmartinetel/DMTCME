gerig_zarate_test <- function(X, y, tabla_critica) {
  bor
  # --- Centrar datos ---
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- y - mean(y)
  
  n <- nrow(Xc)
  k <- ncol(Xc)
  
  # --- SVD ---
  svd_X <- svd(Xc, nu = k, nv = k)
  U <- svd_X$u
  
  # --- Estimar sigma^2 ---
  P_res <- diag(n) - U %*% t(U)
  sigma2_hat <- drop(t(Yc) %*% P_res %*% Yc) / (n - k - 1)
  
  # --- Tabla resultados ---
  resultados_test <- tibble(
    l = integer(),
    s = integer(),
    T_stat = numeric(),
    T_alpha = numeric(),
    decision = character()
  )
  
  ultimo_l_aceptado <- 0
  
  # --- Prueba secuencial ---
  for (l in 1:(k - 1)) {
    
    U_k_s <- U[, (k - l + 1):k, drop = FALSE]
    
    T_stat <- drop(
      t(Yc) %*% U_k_s %*% t(U_k_s) %*% Yc
    ) / (l * sigma2_hat)
    
    T_alpha <- tabla_critica %>%
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
  
  s_sel <- k - ultimo_l_aceptado
  
  list(
    resultados = resultados_test,
    s_seleccionado = s_sel,
    sigma2_hat = sigma2_hat,
    U = U
  )
}