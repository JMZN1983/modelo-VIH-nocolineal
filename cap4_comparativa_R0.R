# ========== Determinista vs Estocástico (R0) — mismos parámetros por escenario ==========
# Se deben cargar los siguientes paquetes
pkgs <- c("ggplot2","patchwork","dplyr","scales","rlang")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(ggplot2); library(patchwork); library(dplyr); library(scales); library(rlang)
set.seed(2025)

# Este segmento de código lo usamos para el panel
pal <- c(negro="#000000", azul="#0072B2", naranja="#E69F00", verde="#009E73",
         cereza="#D55E00", purp="#CC79A7", cielo="#56B4E9", amar="#F0E442", gris="#999999")
tema_q1 <- function(){
  theme_minimal(base_size=13) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=18),
          plot.subtitle=element_text(size=11, margin=margin(b=6)),
          legend.position="bottom", legend.box="horizontal",
          plot.margin=margin(t=8,r=8,b=8,l=8))
}

#Se definen los números reproductivos y funciones de deriva/difusión
R0_RS    <- function(p,B,C,beta,alpha)               (p*B + beta*C)/(alpha + B)                                    # Roberts–Saha (det.)
R0s_Ding <- function(p,B,C,beta0,alpha,gamma0,rho)   ((p-1)*B + (beta0*C - alpha))/gamma0 - (rho^2*C^2)/(2*gamma0) # Ding–Xu–Hu
R0s_NoC  <- function(p,B,C,beta0,alpha,gamma0,rho,s) ((p-1)*B + (beta0*C - alpha))/gamma0 - ((rho^2*C^2+s^2)/(2*gamma0)) # No colineal

F_drift <- function(z,p,B,C,beta0,alpha,gamma0, ...){
  A <- (p-1)*B - gamma0 + (beta0*C - alpha); D <- (beta0*C - alpha)
  A*z - D*z^2
}
G_ding  <- function(z,C,rho, ...)        rho*C*(1-z)*z
G_nocol <- function(z,C,rho,sigma, ...)  { a <- (rho^2)*(C^2); b <- sigma^2; z*sqrt(a*(1-z)^2 + b) }

#  Euler–Maruyama y Euler determinista
em_sim <- function(z0,Tmax,dt,F_fun,G_fun,pars,K=1){
  n <- as.integer(Tmax/dt); t <- seq(0,Tmax,by=dt)
  out <- matrix(NA_real_, n+1, K); out[1,] <- z0
  for(k in 1:K){
    z <- z0
    for(i in 1:n){
      dW <- rnorm(1,0,sqrt(dt))
      z  <- z + rlang::exec(F_fun, z, !!!pars)*dt + rlang::exec(G_fun, z, !!!pars)*dW
      z  <- pmax(0, pmin(1, z))
      out[i+1,k] <- z
    }
  }
  attr(out,"time") <- t; out
}
ode_RS <- function(z0,Tmax,dt,p,B,C,beta,alpha){
  n <- as.integer(Tmax/dt); t <- seq(0,Tmax,by=dt); z <- numeric(n+1); z[1] <- z0
  for(i in 1:n){
    dz <- ( (p-1)*B*z[i] + (beta*C - alpha)*(1-z[i])*z[i] )*dt
    z[i+1] <- pmax(0, pmin(1, z[i] + dz))
  }
  m <- matrix(z, ncol=1); attr(m,"time") <- t; m
}

#  Escalas y leyenda
model_levels <- c("RS","Ding","No-col")
model_cols   <- c("RS"=pal["azul"], "Ding"=pal["naranja"], "No-col"=pal["verde"])
model_lty    <- c("RS"="solid",     "Ding"="dashed",       "No-col"="dotdash")

# Parámetros comunes en los tres modelos dentro de cada escenario
base <- list(p=.20, B=.30, C=1, alpha=.10, beta0=1.20, gamma0=.40, z0=.10, Tmax=70, dt=.1, Nmc=250)

# Escenarios: solo cambian (rho, sigma) para comparar extinción/persistencia entre modelos
S1 <- modifyList(base, list(rho=1.00, sigma=0.60))  # RS>1; Ding<=1; No-col<=1
S2 <- modifyList(base, list(rho=0.60, sigma=0.80))  # RS>1; Ding>1;  No-col<=1
S3 <- modifyList(base, list(rho=0.30, sigma=0.30))  # RS>1; Ding>1;  No-col>1

fmt <- function(x) sprintf("%.2f", x)

# --- Panel: tres modelos (líneas) + bandas IQR para los dos estocásticos
make_panel <- function(sc, tag){
  # RS determinista
  Z_RS <- ode_RS(sc$z0, sc$Tmax, sc$dt, sc$p, sc$B, sc$C, sc$beta0, sc$alpha)
  t    <- attr(Z_RS,"time")
  # Estocásticos (mediana e IQR)
  parsD <- list(p=sc$p,B=sc$B,C=sc$C,beta0=sc$beta0,alpha=sc$alpha,gamma0=sc$gamma0,rho=sc$rho)
  parsN <- list(p=sc$p,B=sc$B,C=sc$C,beta0=sc$beta0,alpha=sc$alpha,gamma0=sc$gamma0,rho=sc$rho,sigma=sc$sigma)
  M_D  <- em_sim(sc$z0, sc$Tmax, sc$dt, F_drift, G_ding,  parsD, K=sc$Nmc)
  M_NC <- em_sim(sc$z0, sc$Tmax, sc$dt, F_drift, G_nocol, parsN, K=sc$Nmc)
  med_D  <- apply(M_D,1,median); q25_D <- apply(M_D,1,quantile,.25); q75_D <- apply(M_D,1,quantile,.75)
  med_NC <- apply(M_NC,1,median); q25_NC<- apply(M_NC,1,quantile,.25); q75_NC<- apply(M_NC,1,quantile,.75)
  
  # R0’s 
  rRS <- R0_RS(sc$p,sc$B,sc$C,sc$beta0,sc$alpha)
  rD  <- R0s_Ding(sc$p,sc$B,sc$C,sc$beta0,sc$alpha,sc$gamma0,sc$rho)
  rNC <- R0s_NoC(sc$p,sc$B,sc$C,sc$beta0,sc$alpha,sc$gamma0,sc$rho,sc$sigma)
  
  # Datos 
  df <- rbind(
    data.frame(time=t, z=as.vector(Z_RS[,1]), modelo="RS"),
    data.frame(time=t, z=med_D,                modelo="Ding"),
    data.frame(time=t, z=med_NC,               modelo="No-col")
  )
  df$modelo <- factor(df$modelo, levels=model_levels)
  
  # Etiqueta por panel
  lbl <- sprintf("RS%s1;  Ding%s1;  No-col%s1\nR0_RS=%s,  R0s_Ding=%s,  R0s=%s",
                 ifelse(rRS>1,">","≤"), ifelse(rD>1,">","≤"), ifelse(rNC>1,">","≤"),
                 fmt(rRS), fmt(rD), fmt(rNC))
  
  ggplot(df, aes(time, z, color=modelo, linetype=modelo)) +
    geom_ribbon(aes(x=time, ymin=q25_D, ymax=q75_D),
                data=data.frame(time=t, q25_D=q25_D, q75_D=q75_D),
                inherit.aes=FALSE, fill=scales::alpha(model_cols["Ding"], .18)) +
    geom_ribbon(aes(x=time, ymin=q25_NC, ymax=q75_NC),
                data=data.frame(time=t, q25_NC=q25_NC, q75_NC=q75_NC),
                inherit.aes=FALSE, fill=scales::alpha(model_cols["No-col"], .14)) +
    geom_line(linewidth=1.2) +
    coord_cartesian(xlim=c(0, sc$Tmax), ylim=c(0,1)) +
    scale_y_continuous(expand=expansion(mult=c(0.02,0.06))) +
    annotate("label", x=0.02*sc$Tmax, y=.96, label=lbl, hjust=0, vjust=1,
             size=3.4, label.size=0, fill=scales::alpha("white",0.75), color="grey20") +
    labs(title=tag, x="Tiempo", y="z(t)") +
    tema_q1()
}

# --- Construcción de paneles y figura final 
p1 <- make_panel(S1,"S1")
p2 <- make_panel(S2,"S2")
p3 <- make_panel(S3,"S3")

fig <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = "Trayectorias — mismos parámetros por escenario (determinista vs estocástico R0)",
    subtitle= sprintf("Parámetros comunes: p=%.2f, B=%.2f, C=%d, γ0=%.2f; α=%.2f; β=%.2f; z0=%.2f",
                      base$p, base$B, base$C, base$gamma0, base$alpha, base$beta0, base$z0)
  ) &
  scale_color_manual(values = setNames(unname(model_cols), model_levels),
                     limits = model_levels, drop = FALSE, name = "Modelo") &
  scale_linetype_manual(values = setNames(unname(model_lty), model_levels),
                        limits = model_levels, drop = FALSE, name = "Modelo") &
  theme(legend.position = "bottom")

# Mostrar en pantalla
fig
