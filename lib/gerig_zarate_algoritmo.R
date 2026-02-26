fit_dtmce_mansfield <- function(X, y,  tabla_critica = NULL, alpha = 0.05 ) {
  # browser()
  X <- as.matrix(X)
  
  # ============================================================
  # 1. CENTRADO
  # ============================================================
  
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- y - mean(y)
  
  n <- nrow(Xc)
  k <- ncol(Xc)
  
  # ============================================================
  # 2. SVD
  # ============================================================
  
  svd_X <- svd(Xc, nu = k, nv = k)
  U <- svd_X$u
  V <- svd_X$v
  d <- svd_X$d
  
  # ============================================================
  # 3. ESTIMACIÓN DE SIGMA
  # ============================================================
  
  P_res <- diag(n) - U %*% t(U)
  sigma2_hat <- drop(t(Yc) %*% P_res %*% Yc) / (n - k - 1)
  
  # ============================================================
  # 4. TORO–WALLACE
  # ============================================================
  
  ultimo_l_aceptado <- 0
  
  for (l in 1:(k - 1)) {
    # browser()
    # 
    U_k_s <- U[, (k - l + 1):k, drop = FALSE]
    
    T_stat <- drop(
      t(Yc) %*% U_k_s %*% t(U_k_s) %*% Yc
    ) / (l * sigma2_hat)
    
    # ============================================================
    # Valor crítico
    # ============================================================
    # browser()
    if (!is.null(tabla_critica)) {
      
      T_alpha <- tabla_critica %>%
        filter(n_k == (n - k - 1)) %>%
        pull(!!sym(paste0("x", l)))
      
    } else {
      
      T_alpha <- qf(1 - alpha, df1 = l, df2 = n - k - 1,ncp = 1)
      
    }
    
    # ============================================================
    
    if (T_stat < T_alpha) {
      ultimo_l_aceptado <- l
    }
  }
  
  s_sel <- max(k- ultimo_l_aceptado, 2)
  
  # ============================================================
  # 5. DTMCE
  # ============================================================
  # browser()
  V_s <- V[, 1:s_sel, drop = FALSE]
  D_s <- diag(d[1:s_sel])
  D_s_inv <- solve(D_s)
  
  X_reducido <- Xc %*% V_s
  modelo_dtmce <- lm(y ~ X_reducido)
  
  gamma_hat <- coef(modelo_dtmce)[-1]
  
  beta_dtmce <- V_s %*% gamma_hat
  
  # ============================================================
  # 6. MANSFIELD
  # ============================================================
  
  calcular_ur <- function(V_s, gamma_hat, D_s_inv, variable_names) {
    
    combs <- combn(seq_len(nrow(V_s)), 1)
    
    ur <- sapply(seq_len(ncol(combs)), function(i) {
      idx <- combs[, i]
      V_r <- V_s[idx, , drop = FALSE]
      
      drop(
        t(gamma_hat) %*% t(V_r) %*%
          (V_r %*% D_s_inv %*% t(V_r)) %*%
          V_r %*% gamma_hat
      )
    })
    
    tibble(
      pair = variable_names,
      ur_value = ur
    ) %>%
      arrange(ur_value)
  }
  
  vars <- colnames(Xc)
  temp_V_s <- V_s
  order_elimin <- character()
  target_vars = s_sel
  while (length(vars) > target_vars) {
    
    ur_table <- calcular_ur(temp_V_s, gamma_hat, D_s_inv, vars)
    
    var_elim <- ur_table$pair[1]
    
    order_elimin <- c(order_elimin, var_elim)
    
    idx_remove <- match(var_elim, vars)
    
    temp_V_s <- temp_V_s[-idx_remove, , drop = FALSE]
    vars <- vars[-idx_remove]
  }
  
  vars_final <- vars
  
  X_final <- Xc[, vars_final, drop = FALSE]
  
  modelo_final <- lm(y ~ ., data = as.data.frame(X_final))
  
  # ============================================================
  # OUTPUT
  # ============================================================
  
  list(
    s_seleccionado = s_sel,
    beta_dtmce = beta_dtmce,
    modelo_dtmce = modelo_dtmce,
    orden_eliminacion = order_elimin,
    variables_finales = vars_final,
    modelo_mansfield = modelo_final
  )
}