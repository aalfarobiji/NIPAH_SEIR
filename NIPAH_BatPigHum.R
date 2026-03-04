
###############################################################################

setwd("D:/NIPAH")

#install.packages("EpiModel")

library("EpiModel")
set.seed(12345)

Nipah <- function(t, t0, parms) { with(as.list(c(t0, parms)), {
  
  # ======================
  # Population totals
  # ======================
  
  num.b <- s_b.num + e_b.num + i_b.num + r_b.num
  num.p <- s_p.num + e_p.num + i_p.num + r_p.num
  num.h <- s_h.num + e_h.num + i_h.num + r_h.num
  
  # ======================
  # Force of Infection
  # ======================
  
  # Bat
  lambda.bb <- beta_BB * (i_b.num/num.b)
  
  # Pig
  lambda.bp <- beta_BP * (i_b.num/num.b)
  lambda.pp <- beta_PP * (i_p.num/num.p)
  
  # Human
  lambda.ph <- beta_PH * (i_p.num/num.p)
  lambda.hh <- beta_HH * (i_h.num/num.h)
  
  # ======================
  # BATS
  # ======================
  dS_b <- -lambda.bb*s_b.num
  dE_b <-  lambda.bb*s_b.num - theta*e_b.num
  dI_b <-  theta*e_b.num - gamma_b*i_b.num
  dR_b <-  gamma_b*(1-delta_b)*i_b.num
  dD_b <-  gamma_b*delta_b*i_b.num

  # ======================
  # PIGS (AMPLIFIER)
  # ======================
  dS_p <- -(lambda.bp + lambda.pp)*s_p.num
  dE_p <-  (lambda.bp + lambda.pp)*s_p.num - omega*e_p.num
  dI_p <-  omega*e_p.num - gamma_p*i_p.num
  dR_p <-  gamma_p*(1-delta_p)*i_p.num
  dD_p <-  gamma_p*delta_p*i_p.num
  
  # ======================
  # HUMANS
  # ======================
  dS_h <- -(lambda.ph + lambda.hh)*s_h.num
  dE_h <-  (lambda.ph + lambda.hh)*s_h.num - sigma*e_h.num 
  
  dI_h <- mu * sigma*e_h.num - rho*i_h.num
  
  dR_h <- rho*(1-delta_h)*i_h.num + (1-mu)*sigma*e_h.num
  dD_h <- rho*delta_h*i_h.num
  
  # ======================
  # Output
  # ======================
  list(c(dS_b, dE_b, dI_b, dR_b, dD_b,
         dS_p, dE_p, dI_p, dR_p, dD_p,
         dS_h, dE_h, dI_h, dR_h, dD_h),
         
         si_p_BP.flow = lambda.bp*s_p.num,
         si_p_PP.flow = lambda.pp*s_p.num,
         si_h_PH.flow = lambda.ph*s_h.num,
         si_h_HH.flow = lambda.hh*s_h.num)
})
}

# ======================
# BASE PARAMETERS
# ======================

param <- param.dcm(
                   # Bat
                   theta   =1/6,
                   gamma_b =1/15,
                   delta_b =0,                     # Kelelawar atau seq(0, 0.1, 0.2)
                   beta_BB =0.01436117,
                   
                   # Pig
                   omega   =1/7,
                   gamma_p =1/10,
                   delta_p =0.7,
                   beta_BP =0.2,
                   beta_PP =0.4,
                   beta_PH =0.4,
  
                   # Human
                   sigma   =1/9,
                   rho     =1/6,
                   mu      =0.3,                  # c(0.1, 0.3, 0.5), #seq(0.0, 0.3, 0.03),  # Vary Mu from 0 to 0.3 in increments of 0.03  # Kerala 0.076
                   #delta1  =1,
                   delta_h =0.71315524,
                   beta_HH =0.69733082)
                   #beta2_HH=0.04760139)

init <- init.dcm(
# Bat
s_b.num = 500, e_b.num=0, i_b.num = 0, r_b.num = 0, d_b.num = 0,
# Pig
s_p.num = 500, e_p.num=0, i_p.num = 0, r_p.num = 0, d_p.num = 0,
# Human
s_h.num = 999, e_h.num=0, i_h.num = 1, r_h.num = 0, d_h.num = 0)

control <- control.dcm(nsteps = 365, new.mod = Nipah)

# ======================
# BASELINE SIMULATION
# ======================

sim <- dcm(param, init, control)

sim <- mutate_epi(sim,
                  # Populations Total
                  b.total = s_b.num + e_b.num + i_b.num + r_b.num + d_b.num,
                  p.total = s_p.num + e_p.num + i_p.num + r_p.num + d_p.num,
                  h.total = s_h.num + e_h.num + i_h.num + r_h.num + d_h.num,

                  i_b.total = i_b.num,
                  i_p.total = i_p.num,
                  i_h.total = i_h.num,

                  # Cumulative Infection Flows
                  cum_BP = cumsum(si_p_BP.flow),
                  cum_PP = cumsum(si_p_PP.flow),
                  cum_PH = cumsum(si_h_PH.flow),
                  cum_HH = cumsum(si_h_HH.flow)) 

sim <- mutate_epi(sim,
                  cum_total_h = cum_PH + cum_HH)          # Total Outbreak pada Manusia

#sim <- mutate_epi(sim,                                                          
                  #prop_spillover = ifelse(cum_total_h>0,
                                          #cum_PH/cum_total_h,0)) # Proporsi Outbreak dari Hewan

sim <- mutate_epi(sim,                                                          
                  cum_death_h = d_h.num)

sim <- mutate_epi(sim,
                  cum_infec_h = init$s_h.num - s_h.num)

c(
  cumulative_infection = tail(sim$epi$cum_infec_h,1),
  cumulative_death     = tail(sim$epi$cum_death_h,1))

#sum(sim$epi$cum_death_h)
#sum(sim$epi$i_h.total)
#tail(sim$epi$cum_death_h)
#tail(sim$epi$prop_spillover) # Melihat hasil hari terakhir
#summary(sim$epi$prop_spillover)
#tail(sim$epi$i_h.total)

# ======================
# R0 SECTION (Integrated)
# ======================

# Ambil parameter dari object param
with(as.list(param), {
  
  # ----------------------
  # R0 BAT
  # ----------------------
  R0_bat <<- beta_BB / gamma_b
  
  # ----------------------
  # R0 PIG (within pig)
  # ----------------------
  R0_pig <<- beta_PP / gamma_p
  
  # ----------------------
  # R0 HUMAN (HH only)
  # ----------------------
  R0_human <<- (beta_HH * mu) / rho
  
  # ----------------------
  # FULL MULTI-HOST R0 via NGM
  # ----------------------
  
  # F matrix (new infections)
  F_mat <- matrix(c(
    beta_BB, 0,        0,
    beta_BP, beta_PP,  0,
    0,       beta_PH,  beta_HH * mu
  ), nrow = 3, byrow = TRUE)
  
  # V matrix (transition)
  V_mat <- matrix(c(
    gamma_b, 0,       0,
    0,       gamma_p, 0,
    0,       0,       rho
  ), nrow = 3, byrow = TRUE)
  
  # Next Generation Matrix
  K_mat <- F_mat %*% solve(V_mat)
  
  eigs <- eigen(K_mat)$values
  
  R0_full <<- max(Re(eigs))
  
  dominant_driver <<- c("Bat", "Pig", "Human")[which.max(Re(eigs))]
})

# ======================
# PRINT RESULTS
# ======================

cat("\n==============================\n")
cat("R0 Bat                =", round(R0_bat,3), "\n")
cat("R0 Pig (within pig)   =", round(R0_pig,3), "\n")
cat("R0 Human (HH only)    =", round(R0_human,3), "\n")
cat("Full Multi-Host R0    =", round(R0_full,3), "\n")
cat("Dominant Driver       =", dominant_driver, "\n")
cat("==============================\n\n")

# ======================
# Perbandingan Angka Infeksi
# ======================
par(bty = "l", mar = c(5,5,4,2))       

time <- 1:365

data_mat <- cbind(sim$epi$i_b.total,
                  sim$epi$i_p.total,
                  sim$epi$i_h.total)

matplot(time,
        data_mat,
        type = "l",
        lty = 1,
        lwd = 2,
        col = c("blue","red","darkgreen"),
        xlab = "Time (Days)",
        ylab = "Number",
        main = "Comparison of Infections",
        ylim = c(0,5))

cols <- c("blue","red","darkgreen")

for(i in 1:ncol(data_mat)){
  peak_val <- max(data_mat[,i])
  peak_day <- which.max(data_mat[,i])

  abline(v=peak_day,
         col=adjustcolor(cols[i], alpha.f=0.6),
         lty=2,
         lwd=1)
  
  abline(h=peak_val,
         col=adjustcolor(cols[i], alpha.f=0.6),
         lty=2,
         lwd=1)
  
  text(peak_day,
       peak_val,
       labels=paste0("(", peak_day, ", ", round(peak_val,1), ")"),
       pos=3,
       col=cols[i],
       cex=0.9)
}

legend("topright",
       inset = 0.03,   
       bty = "n",      
       legend = c("Bat","Pig","Human"),
       col = c("blue","red","darkgreen"),
       lty = 1,
       lwd = 2)

# ======================
# Perbandingan Angka Kematian
# ======================
par(bty = "l", mar = c(5,5,4,2))       

time <- 1:365

data_mat <- cbind(sim$epi$d_b.num,
                  sim$epi$d_p.num,
                  sim$epi$d_h.num)

matplot(time,
        data_mat,
        type = "l",
        lty = 1,
        lwd = 2,
        col = c("blue","red","darkgreen"),
        xlab = "Time (Days)",
        ylab = "Number",
        main = "Comparison of Death",
        ylim = c(0,10))

cols <- c("blue","red","darkgreen")

for(i in 1:ncol(data_mat)){
  peak_val <- max(data_mat[,i])
  peak_day <- which.max(data_mat[,i])
  
  abline(v=peak_day,
         col=adjustcolor(cols[i], alpha.f=0.6),
         lty=2,
         lwd=1)
  
  abline(h=peak_val,
         col=adjustcolor(cols[i], alpha.f=0.6),
         lty=2,
         lwd=1)
  
  text(peak_day,
       peak_val,
       labels=paste0("(", peak_day, ", ", round(peak_val,1), ")"),
       pos=3,
       col=cols[i],
       cex=0.9)
}

legend("topright",
       inset = 0.03,   
       bty = "n",      
       legend = c("Bat","Pig","Human"),
       col = c("blue","red","darkgreen"),
       lty = 1,
       lwd = 2)

# ======================
# BOXPLOT R0 COMPARISON
# ======================

R0_values <- data.frame(
  Bat    = R0_bat,
  Pig    = R0_pig,
  Human  = R0_human,
  Full   = R0_full
)

par(mar=c(5,6,4,2), bty="l")

boxplot(R0_values,
        col=c("blue","red","darkgreen","purple"),
        ylab=expression(R[0]),
        main=expression("Comparison of Basic Reproduction Numbers ("*R[0]*")"))

# Garis threshold R0 = 1
abline(h=1, col="black", lty=2, lwd=2)

# Tambahkan nilai numerik
points(1:4,
       c(R0_bat, R0_pig, R0_human, R0_full),
       pch=19)

# ======================
# R0 ANALYSIS 
# ======================

# Sequence mu
mu_seq <- seq(0,0.5,by=0.001)

# Formal R0 (NGM-based)
R0_seq <- (param$beta_HH * mu_seq) / param$rho

# Critical mu when R0 = 1
mu_critical <- param$rho / param$beta_HH

# ======================
# PLOT R0 vs MU
# ======================

plot(mu_seq, R0_seq,
     type="l", lwd=2, col="purple",
     main=expression("Basic Reproduction Number ("*R[0]*") vs "*mu),
     xlab=expression(mu),
     ylab=expression(R[0]),
     ylim=c(0, max(R0_seq)))

# Threshold R0 = 1
abline(h=1, col="red", lty=2, lwd=2)

# Critical mu
abline(v=mu_critical, col="red", lty=2, lwd=2)

# Kerala 2018 reference
mu_kerala <- 0.07621085
R0_kerala <- (param$beta_HH * mu_kerala) / param$rho

abline(v=mu_kerala, col="blue", lwd=2, lty=2)
abline(h=R0_kerala, col="blue", lwd=2, lty=2)

# Baseline (mu = 0.3)
mu_baseline <- param$mu
R0_baseline <- (param$beta_HH * mu_baseline) / param$rho

abline(v=mu_baseline, col="darkgreen", lwd=2)
abline(h=R0_baseline, col="darkgreen", lwd=2)

legend(x = 0.0, y = 2.0,
       legend=c(expression(R[0]),
                "R0 = 1",
                "Kerala 2018",
                "Baseline (mu=0.3)"),
       col=c("purple","red","blue","darkgreen"),
       lwd=2, lty=c(1,2,2,1),
       bty="0",
       bg  = "white",
       cex = 0.8)

############################################################################