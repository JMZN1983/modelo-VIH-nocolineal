# =================== VIH con dos ruidos no colineales — Extinción, Persistencia, Umbral ===================
# Paquetes necesarios
pkgs <- c("ggplot2","patchwork","dplyr","knitr","scales","rlang")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(ggplot2); library(patchwork); library(dplyr); library(knitr); library(scales); library(rlang)
set.seed(2025)

# --------------------------- Estilo ---------------------------
pal <- c(negro="#000000", azul="#0072B2", naranja="#E69F00", verde="#009E73",
         cereza="#D55E00", purp="#CC79A7", cielo="#56B4E9", amar="#F0E442", gris="#999999")
tema_q1 <- function(){
  theme_minimal(base_size=13) +
    theme(panel.grid.minor=element_blank(),
          legend.position="top",
          strip.text=element_text(face="bold"),
          plot.title=element_text(face="bold", size=16))
}

# --------------------------- Número reproductivo estocástico (no colineal) ---------------------------
R0s_NoC <- function(p,B,C,beta0,alpha,gamma0,rho,sigma){
  ((p-1)*B + (beta0*C - alpha))/gamma0 - ((rho^2*C^2 + sigma^2)/(2*gamma0))
}

# --------------------------- Deriva y difusión del modelo con dos ruidos no colineales ---------------------------
F_drift <- function(z,p,B,C,beta0,alpha,gamma0, ...){
  A <- (p-1)*B - gamma0 + (beta0*C - alpha)   # término lineal
  D <- (beta0*C - alpha)                      # término cuadrático (logístico)
  A*z - D*z^2
}
G_nocol <- function(z,C,rho,sigma, ...){
  a <- (rho^2)*(C^2); b <- sigma^2
  z*sqrt(a*(1-z)^2 + b)
}

# --------------------------- Integrador Euler–Maruyama ---------------------------
em_sim <- function(z0,Tmax,dt,F_fun,G_fun,pars,K=1){
  n <- as.integer(Tmax/dt); t <- seq(0,Tmax,by=dt)
  out <- matrix(NA_real_, n+1, K); out[1,] <- z0
  for(k in 1:K){
    z <- z0
    for(i in 1:n){
      dW <- rnorm(1,0,sqrt(dt))
      z  <- z + rlang::exec(F_fun, z, !!!pars)*dt + rlang::exec(G_fun, z, !!!pars)*dW
      z  <- max(0, min(1, z))  # reflexión blanda en [0,1]
      out[i+1,k] <- z
    }
  }
  attr(out,"time") <- t; out
}

# --------------------------- Densidad estacionaria π(z) (no colineal) ---------------------------
pi_explicit_core <- function(z,p,B,C,beta0,alpha,gamma0,rho,sigma){
  a <- rho^2*C^2; b <- sigma^2; k2 <- a + b
  A <- (p-1)*B - gamma0 + (beta0*C - alpha); D <- (beta0*C - alpha)
  H <- sqrt(a*(1-z)^2 + b); eta <- if(b==0) Inf else sqrt(a/b)
  base <- z^( (2*A/k2) - 2 ) * H^( -(2*A/k2 + 2) )
  if(is.finite(eta)){
    base * exp( -(2/sqrt(a*b))*(D + (a*A/k2))*atan( eta*(1-z) ) )
  } else NA_real_
}
pi_norm <- function(p,B,C,beta0,alpha,gamma0,rho,sigma){
  f <- function(x) pi_explicit_core(x,p,B,C,beta0,alpha,gamma0,rho,sigma)
  val <- integrate(function(x) pmax(1e-300, f(x)), lower=1e-6, upper=1-1e-6,
                   rel.tol=1e-7, subdivisions=5000)$value
  1/val
}
pi_eval <- function(z,p,B,C,beta0,alpha,gamma0,rho,sigma){
  K <- pi_norm(p,B,C,beta0,alpha,gamma0,rho,sigma)
  K * pi_explicit_core(z,p,B,C,beta0,alpha,gamma0,rho,sigma)
}

# --------------------------- Gráficos auxiliares ---------------------------
g_trayectorias <- function(caso, T_max, dt=.1, N_mc=220, N_plot=60, mediana=TRUE, IQR=TRUE, salto_x=10){
  pars <- list(p=caso$p,B=caso$B,C=caso$C,beta0=caso$beta0,alpha=caso$alpha,gamma0=caso$gamma0,
               rho=caso$rho,sigma=caso$sigma)
  M <- em_sim(caso$z0, T_max, dt, F_drift, G_nocol, pars, K=N_mc)
  t <- attr(M,"time"); cada <- 10; idx <- seq(1,length(t),by=cada)
  keep <- sample(seq_len(ncol(M)), min(N_plot, ncol(M))); Mv <- M[,keep,drop=FALSE]
  dfL <- data.frame(tiempo=rep(t[idx], times=ncol(Mv)),
                    z=as.vector(Mv[idx,,drop=FALSE]),
                    camino=rep(seq_len(ncol(Mv)), each=length(idx)))
  zt <- split(as.vector(M[idx,,drop=FALSE]), rep(seq_along(idx), each=ncol(M)))
  med <- sapply(zt, median); q25 <- sapply(zt, quantile, 0.25); q75 <- sapply(zt, quantile, 0.75)
  dfS <- data.frame(tiempo=t[idx], med=med, q25=q25, q75=q75)
  ggplot() +
    {if(IQR) geom_ribbon(data=dfS, aes(x=tiempo,ymin=q25,ymax=q75), fill="grey85")} +
    geom_line(data=dfL, aes(tiempo,z,group=camino), color=pal["gris"], alpha=.55, linewidth=.35) +
    {if(mediana) geom_line(data=dfS, aes(tiempo,med), color=pal["azul"], linewidth=1.05)} +
    scale_x_continuous(limits=c(0,T_max), breaks=seq(0,T_max,by=salto_x)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,by=.2)) +
    labs(x="Tiempo", y="z(t)") + tema_q1()
}
g_hist_vs_pi <- function(caso, T_max=300, dt=.1, N_mc=160, frac_cola=.5, salto_cola=20){
  if(caso$sigma<=0 || R0s_NoC(caso$p,caso$B,caso$C,caso$beta0,caso$alpha,caso$gamma0,caso$rho,caso$sigma) <= 1){
    return(ggplot() + annotate("text", x=.5,y=.5, label="Sin densidad estacionaria (R0s ≤ 1 o σ=0)", size=5) +
             labs(x="z",y="Densidad") + tema_q1())
  }
  pars <- list(p=caso$p,B=caso$B,C=caso$C,beta0=caso$beta0,alpha=caso$alpha,gamma0=caso$gamma0,
               rho=caso$rho,sigma=caso$sigma)
  M <- em_sim(caso$z0, T_max, dt, F_drift, G_nocol, pars, K=N_mc)
  n <- nrow(M); i0 <- max(2, floor(n*(1-frac_cola))); idx <- seq(i0, n, by=salto_cola)
  zc <- as.numeric(M[idx,,drop=FALSE])
  zgrid <- seq(0.001,0.999,length.out=600)
  dpi <- data.frame(z=zgrid, pi=pi_eval(zgrid, caso$p,caso$B,caso$C,caso$beta0,caso$alpha,caso$gamma0,caso$rho,caso$sigma))
  mu  <- integrate(function(x) x*pi_eval(x,caso$p,caso$B,caso$C,caso$beta0,caso$alpha,caso$gamma0,caso$rho,caso$sigma),
                   0,1, rel.tol=1e-6, subdivisions=4000)$value
  m2  <- integrate(function(x) x^2*pi_eval(x,caso$p,caso$B,caso$C,caso$beta0,caso$alpha,caso$gamma0,caso$rho,caso$sigma),
                   0,1, rel.tol=1e-6, subdivisions=4000)$value
  sdT <- sqrt(max(0, m2 - mu^2)); muE <- mean(zc)
  ggplot(data.frame(z=zc), aes(z)) +
    geom_histogram(aes(y=after_stat(..density..)), bins=60, fill=pal["amar"], color="white") +
    geom_line(data=dpi, aes(z,pi), color=pal["azul"], linewidth=1.05) +
    geom_vline(xintercept=mu, color=pal["cereza"], linewidth=1.05) +
    geom_vline(xintercept=c(mu-sdT, mu+sdT), color=pal["cereza"], linetype="dotted", linewidth=.9) +
    geom_vline(xintercept=muE, color=pal["negro"], linetype="longdash", linewidth=.9) +
    coord_cartesian(xlim=c(0,1)) + labs(x="z", y="Densidad") + tema_q1()
}

# --------------------------- Casos (parámetros) ---------------------------
mk_caso <- function(nombre,p=.20,B=.30,C=1,alpha=.10,beta0=1.20,gamma0=.50,rho=.40,sigma=.40,z0=.10){
  list(nombre=nombre,p=p,B=B,C=C,alpha=alpha,beta0=beta0,gamma0=gamma0,rho=rho,sigma=sigma,z0=z0)
}
# A: Extinción
casos_A <- list(
  mk_caso("A1", beta0=.60,  p=.20,B=.55, alpha=.40, gamma0=.45, rho=.60, sigma=.40, z0=.10),
  mk_caso("A2", beta0=.55,  p=.20,B=.55, alpha=.45, gamma0=.45, rho=.70, sigma=.40, z0=.20),
  mk_caso("A3", beta0=.65,  p=.20,B=.60, alpha=.45, gamma0=.45, rho=.70, sigma=.60, z0=.80),
  mk_caso("A4", beta0=.70,  p=.20,B=.60, alpha=.45, gamma0=.45, rho=.80, sigma=.30, z0=.50),
  mk_caso("A5", beta0=.50,  p=.20,B=.50, alpha=.45, gamma0=.45, rho=.85, sigma=.50, z0=.30),
  mk_caso("A6", beta0=.68,  p=.20,B=.60, alpha=.45, gamma0=.45, rho=.85, sigma=.40, z0=.70)
)
# B: Persistencia (forzamos algunos a R0s_NoC >= 1.15)
casos_B_raw <- list(
  mk_caso("B1", beta0=1.20, p=.20,B=.30, alpha=.10, gamma0=.50, rho=.35, sigma=.40, z0=.10),
  mk_caso("B2", beta0=1.50, p=.20,B=.30, alpha=.10, gamma0=.50, rho=.40, sigma=.40, z0=.30),
  mk_caso("B3", beta0=1.20, p=.20,B=.20, alpha=.10, gamma0=.40, rho=.40, sigma=.40, z0=.50),
  mk_caso("B4", beta0=1.20, p=.20,B=.30, alpha=.10, gamma0=.50, rho=.28, sigma=.45, z0=.70),
  mk_caso("B5", beta0=1.40, p=.20,B=.30, alpha=.10, gamma0=.50, rho=.40, sigma=.40, z0=.90),
  mk_caso("B6", beta0=1.35, p=.20,B=.25, alpha=.10, gamma0=.45, rho=.40, sigma=.60, z0=.20)
)
ajusta_persistencia <- function(c, objetivo=1.15){
  R <- R0s_NoC(c$p,c$B,c$C,c$beta0,c$alpha,c$gamma0,c$rho,c$sigma)
  if(is.na(R) || R < objetivo){
    delta <- (objetivo - R)*c$gamma0
    c$beta0 <- c$beta0 + 1.05*delta
  }
  c
}
casos_B <- Map(function(x){ if(x$nombre %in% c("B1","B4","B6")) ajusta_persistencia(x,1.15) else x }, casos_B_raw)

# D: Umbral (R0s_NoC = 1) — despeje de beta*
mk_umbral <- function(p=.20,B=.30,C=1,alpha=.10,gamma0=.50,rho=.50,sigma=.40,z0=.50){
  beta_star <- (gamma0 + alpha - (p-1)*B + (rho^2*C^2 + sigma^2)/2)/C
  mk_caso("Umbral", p=p,B=B,C=C,alpha=alpha,beta0=beta_star,gamma0=gamma0,rho=rho,sigma=sigma,z0=z0)
}
caso_D <- mk_umbral()

# --------------------------- Figuras principales ---------------------------
# Figura 1: Extinción (trayectorias, T=25)
gA <- lapply(casos_A, g_trayectorias, T_max=25, dt=.1, N_mc=220, N_plot=60, mediana=TRUE, IQR=TRUE, salto_x=5)
fig1 <- (gA[[1]]|gA[[2]]|gA[[3]])/(gA[[4]]|gA[[5]]|gA[[6]]) +
  plot_annotation(title="Figura 1. Extinción — 6 casos (trayectorias; Tmáx=25)") & tema_q1()

# Figura 2: Persistencia (arriba trayectorias; abajo densidades estacionarias)
gB  <- lapply(casos_B, g_trayectorias, T_max=300, dt=.1, N_mc=220, N_plot=60, mediana=TRUE, IQR=TRUE, salto_x=100)
gPi <- lapply(casos_B, g_hist_vs_pi, T_max=300, dt=.1, N_mc=160, frac_cola=.5, salto_cola=20)
fig2 <- (gB[[1]]|gB[[2]]|gB[[3]])/(gPi[[1]]|gPi[[2]]|gPi[[3]])/(gB[[4]]|gB[[5]]|gB[[6]])/(gPi[[4]]|gPi[[5]]|gPi[[6]]) +
  plot_annotation(title="Figura 2. Persistencia — trayectorias (arriba) y densidades estacionarias (abajo)") & tema_q1()

# Figura 3: Umbral (R0s = 1)
g3_arriba <- g_trayectorias(caso_D, T_max=300, dt=.1, N_mc=220, N_plot=60, mediana=TRUE, IQR=TRUE, salto_x=100) +
  labs(title="Figura 3 (arriba). Umbral (R0s = 1) — trayectorias; Tmáx=300")
g3_abajo <- {
  parsD <- list(p=caso_D$p,B=caso_D$B,C=caso_D$C,beta0=caso_D$beta0,alpha=caso_D$alpha,gamma0=caso_D$gamma0,
                rho=caso_D$rho,sigma=caso_D$sigma)
  t <- seq(0,300,by=.1); n <- length(t)
  M <- em_sim(caso_D$z0, 300, .1, F_drift, G_nocol, parsD, K=200)
  idx <- seq(max(2, floor(n*.8)), n, by=5); zf <- as.numeric(M[idx,,drop=FALSE])
  ggplot(data.frame(z=zf), aes(z)) +
    geom_histogram(aes(y=after_stat(..density..)), bins=60, fill=pal["amar"], color="white") +
    coord_cartesian(xlim=c(0,1)) +
    labs(title="Figura 3 (abajo). Umbral — densidad empírica final (recurrente nulo, no estacionaria)",
         x="z", y="Densidad") + tema_q1()
}
fig3 <- g3_arriba / g3_abajo

# --------------------------- Tabla LaTeX (parámetros por figura/panel) ---------------------------
to_row <- function(lista, pref,Tmax,dt){
  do.call(rbind, lapply(seq_along(lista), function(i){
    c <- lista[[i]]
    data.frame(Panel=sprintf("%s%d",pref,i),
               p=c$p,B=c$B,C=c$C,alpha=c$alpha,beta0=c$beta0,gamma0=c$gamma0,rho=c$rho,sigma=c$sigma,z0=c$z0,
               R0s_NoCol=round(R0s_NoC(c$p,c$B,c$C,c$beta0,c$alpha,c$gamma0,c$rho,c$sigma),3),
               T_max=Tmax, dt=dt, check.names=FALSE)
  }))
}
tabla_param <- dplyr::bind_rows(
  to_row(casos_A,"1-",25,.1),
  to_row(casos_B,"2-",300,.1),
  to_row(list(caso_D),"3-",300,.1)
)
tabla_latex <- knitr::kable(
  tabla_param, format="latex", booktabs=TRUE,
  caption="Par\\'ametros de simulaci\\'on por figura/panel (modelo VIH con dos ruidos no colineales).",
  label="tab:parametros_VIH",
  col.names=c("Panel","$p$","$B$","$C$","$\\alpha$","$\\beta_0$","$\\gamma_0$","$\\rho$","$\\sigma$",
              "$z_0$","$R_0^{\\mathrm{s}}$","$T_{\\max}$","$\\Delta t$")
)

# --------------------------- Exportación (PNG y PDF) ---------------------------
dir_out <- "D:/Users/Jose/Downloads/VIH"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
dev_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ts <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
guardar <- function(g, nombre, w, h, dpi=300){
  p_pdf <- file.path(dir_out, paste0(nombre,"_",ts(),".pdf"))
  p_png <- file.path(dir_out, paste0(nombre,"_",ts(),".png"))
  try(ggsave(p_pdf, g, width=w, height=h, device=dev_pdf), silent=TRUE)
  ggsave(p_png, g, width=w, height=h, dpi=dpi)
  message("Guardado: ", p_png)
}
guardar(fig1, "Figura_1_Extincion_T25", 11, 7)
guardar(fig2, "Figura_2_Persistencia_Pareada_T300", 11, 14)
guardar(fig3, "Figura_3_Umbral_R0s1_T300", 11, 10)
writeLines(tabla_latex, file.path(dir_out, "parametros_figuras_VIH.tex"))
# ---- FIN ----
