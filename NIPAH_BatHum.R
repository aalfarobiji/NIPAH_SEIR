
#########################################################################

setwd("D:/NIPAH")

#install.packages("EpiModel")

library("EpiModel")
set.seed(12345)

Nipah <- function(t, t0, parms) { with(as.list(c(t0, parms)), {
  
  # ======================
  # Population totals
  # ======================
  
  num.b <- s_b.num + e_b.num + i_b.num + r_b.num
  num.h <- s_h.num + e_h.num + i_h.num + r_h.num
  
  # ======================
  # Force of Infection
  # ======================
  
  lambda.bb <-  beta_BB * (i_b.num/num.b)
  lambda.bh <-  beta_BH * (i_b.num/num.b)
  lambda.hh <-  beta_HH * (i_h.num/num.h)
  
  # ======================
  # Differential equations
  # ======================
  # Bats
  dS_b <- -lambda.bb*s_b.num
  dE_b <- lambda.bb*s_b.num - theta*e_b.num
  dI_b <- theta*e_b.num - gamma*i_b.num
  dR_b <- gamma*(1-delta_b)*i_b.num
  dD_b <- gamma*delta_b*i_b.num
  
  # Human
  dS_h <- -lambda.hh*s_h.num - lambda.bh*s_h.num
  dE_h <- (lambda.hh + lambda.bh)*s_h.num - sigma*e_h.num
  dI_h <- mu*sigma*e_h.num - rho*i_h.num
  dR_h <- rho*(1-delta_h)*i_h.num + (1-mu)*sigma*e_h.num  
  dD_h <- rho*delta_h*i_h.num
  
  # ======================
  # Output
  # ======================
  list(c(dS_b, dE_b, dI_b, dR_b, dD_b,
         dS_h, dE_h, dI_h, dR_h, dD_h,
         si_h_BH.flow = lambda.bh*s_h.num,
         si_h_HH.flow = lambda.hh*s_h.num))
})
}

# ======================
# BASE PARAMETERS
# ======================

param <- param.dcm(theta=1/6,
                   gamma=1/15,
                   sigma=1/9,
                   rho=1/6,
                   mu=c(0.1, 0.3, 0.5),            # seq(0.0, 0.3, 0.03), # Vary Mu from 0 to 0.3 in increments of 0.03
                   #delta1=1,
                   delta_h=0.71315524,
                   delta_b=0,                     # Kelelawar atau c(0, 0.1, 0.2)
                   beta_BB=0.01436117,
                   beta_BH=0.11595493,
                   beta_HH=0.69733082)
                   #beta2_HH=0.04760139)

init <- init.dcm(s_b.num = 999, e_b.num=0, i_b.num = 1, r_b.num = 0, d_b.num = 0,
                 s_h.num = 1000, e_h.num=0, i_h.num = 0, r_h.num = 0, d_h.num = 0,
                 si_h_BH.flow=0, si_h_HH.flow=0)

control <- control.dcm(nsteps = 365, new.mod = Nipah)

# ======================
# BASELINE SIMULATION
# ======================

sim <- dcm(param, init, control)

#sim <- mutate_epi(sim, h.total = s_h.num + e_h.num + i1_h.num + i2_h.num + r_h.num) # Total number of humans
#sim <- mutate_epi(sim, i_h.total = i1_h.num + i2_h.num)  

sim <- mutate_epi(sim,
                  h.total = s_h.num + e_h.num + i_h.num + r_h.num,
                  i_h.total = i_h.num,
                  b.total = s_b.num + e_b.num + i_b.num + r_b.num,
                  i_b.total = i_b.num)

sim <- mutate_epi(sim,
                  cum_BH = cumsum(si_h_BH.flow),
                  cum_HH = cumsum(si_h_HH.flow)) # Cumulative Spillover (1)

sim <- mutate_epi(sim,
                  cum_total = cum_BH + cum_HH) # Total Outbreak pada Manusia

sim <- mutate_epi(sim,
                  prop_spillover = cum_BH / cum_total) # Proporsi Outbreak dari Hewan

#tail(sim$epi$prop_spillover) # Melihat hasil hari terakhir
#summary(sim$epi$prop_spillover)
#summary(sim$epi$i_h.total)

# ======================
# Perbandingan Angka Infeksi
# ======================

par(bty = "l")              
par(mar = c(5,5,4,2))       

time <- 1:365

matplot(time,
        sim$epi$i_h.total,
        type = "l",
        lty = 1,
        lwd = 2,
        col = c("blue","red","darkgreen"),
        xlab = "Time (Days)",
        ylab = "Number",
        main = "Human Infection Under Different Mu",
        ylim = c(0,60))

cols <- c("blue","red","darkgreen")

for(i in 1:ncol(sim$epi$i_h.total)){
  peak_val <- max(sim$epi$i_h.total[,i])
  peak_day <- which.max(sim$epi$i_h.total[,i])
  
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
       legend = c("0.1", "0.3", "0.5"),
       col = c("blue","red","darkgreen"),
       lty = 1,
       lwd = 2)

# ======================
# Perbandingan Angka Kematian
# ======================

par(bty = "l")              
par(mar = c(5,5,4,2))       

time <- 1:365

matplot(time,
        sim$epi$d_h.num,
        type = "l",
        lty = 1,
        lwd = 2,
        col = c("blue","red","darkgreen"),
        xlab = "Time (Days)",
        ylab = "Number",
        main = "Human Death Under Different Mu",
        ylim = c(0,600))

cols <- c("blue","red","darkgreen")

for(i in 1:ncol(sim$epi$d_h.num)){
  peak_val <- max(sim$epi$d_h.num[,i])
  peak_day <- which.max(sim$epi$d_h.num[,i])
  
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
       legend = c("0.1", "0.3", "0.5"),
       col = c("blue","red","darkgreen"),
       lty = 1,
       lwd = 2)

     #plot(sim,
     #y = "i_h.num",
     #main = "SS VS NS",
     #legend = "full",
     #ylim = c(0,60))      

     #plot(sim,
     #y = c("cum_BH", "cum_HH"),
     #main = "Spillover VS HH",
     #legend = "full",
     #ylim = c(0,3))      # Plot Cumulative Spillover 

     #plot(sim,
     #y = "prop_spillover",
     #main = "Proportion of Spillover",
     #legend = "full",
     #ylim = c(0,0.4),
     #xlim = c(0,400))   # Plot Proporsi Spillover

     #plot(sim, y = "i_h.total")

     #write.csv(sim, "NIPAH_summary_results.csv", row.names = FALSE)

     #legend(x = 10, y = 1000,
       #legend = paste("run", 1:11, sep=""),
       #col = 1:11,
       #lty = 1)

# ======================
# RANGE MU
# ======================

mu_seq <- seq(0, 0.6, by = 0.01)

beta_HH_val <- param$beta_HH[1]
rho_val     <- param$rho[1]

# R0 kontinu berdasarkan model
R0_seq <- (beta_HH_val * mu_seq) / rho_val

# 3 skenario mu dari simulasi
mu_values <- param$mu
R0_values <- (beta_HH_val * mu_values) / rho_val

# ======================
# HITUNG THRESHOLD OTOMATIS
# ======================

mu_critical <- rho_val / beta_HH_val   # R0 = 1 condition

# ======================
# KERALA 2018 (FIXED REFERENCE)
# ======================

mu_kerala  <- 0.07621085
R0_kerala  <- 0.60   # nilai literature

# ======================
# PLOT
# ======================

par(bty="l")
par(mar=c(5,5,4,2))

plot(mu_seq, R0_seq,
     type="l",
     lwd=3,
     col="purple",
     main=expression("Basic Reproduction Number ("*R[0]*") vs "*mu),
     xlab=expression(mu),
     ylab=expression(R[0]),
     ylim=c(0, max(R0_seq)+0.5))

# ----------------------
# THRESHOLD (MODEL-BASED)
# ----------------------

abline(h=1, col="red", lwd=2, lty=2)
abline(v=mu_critical, col="red", lwd=2, lty=2)

#text(mu_critical, 1,
     #labels=paste0("mu* = ", round(mu_critical,3)),
     #pos=4,
     #col="red")

# ----------------------
# KERALA 2018 (BIRU)
# ----------------------

abline(h=R0_kerala, col="blue", lwd=2, lty=2)
abline(v=mu_kerala, col="blue", lwd=2, lty=2)

#points(mu_kerala, R0_kerala,
       #pch=17,
       #cex=2,
       #col="blue")

# ----------------------
# 3 SKENARIO SIMULASI
# ----------------------

points(mu_values, R0_values,
       pch=19,
       cex=2,
       col=c("darkgreen","orange","brown"))

#text(mu_values,
     #R0_values,
     #labels=paste0("mu=", mu_values),
     #pos=3)

# ----------------------
# LEGEND
# ----------------------

legend(x = 0.0, y = 3.0,
       legend = c(expression(R[0]~"(model)"),
                  "Threshold (R0 = 1)",
                  "Kerala 2018 (R0 = 0.60)",
                  "Simulation: mu = 0.1",
                  "Simulation: mu = 0.3",
                  "Simulation: mu = 0.5"),
       col = c("purple",
               "red",
               "blue",
               "darkgreen",
               "orange",
               "brown"),
       lwd = c(3, 2, 2, NA, NA, NA),
       lty = c(1, 2, 2, NA, NA, NA),
       pch = c(NA, NA, NA, 19, 19, 19),
       bty = "o",
       bg  = "white",
       cex = 0.8)

#########################################################################

# ======================
# INTERVENTION PARAMETERS
# ======================

# Misalnya:
# - 20% reduction in human-to-human transmission (Isolasi, APD, tracing, hospital IPC)
# - 10% reduction in progression to infectious Early detection → lebih banyak kasus abortive / tidak jadi infectious

beta_HH_int <- param$beta_HH * 0.8
mu_int      <- param$mu * 0.9

param_int <- param.dcm(theta=1/6,
                       gamma=1/15,
                       sigma=1/9,
                       rho=1/6,
                       mu=mu_int,
                       delta_h=0.71315524,
                       delta_b=0,
                       beta_BB=0.01436117,
                       beta_BH=0.11595493,
                       beta_HH=beta_HH_int)

# ======================
# RUN INTERVENTION SIMULATION
# ======================

sim_int <- dcm(param_int, init, control)

sim_int <- mutate_epi(sim_int,
                      i_h.total = i_h.num)

# ======================
# COMPARISON PLOT
# ======================

par(bty="l")
par(mar=c(5,5,4,2))

matplot(time,
        sim$epi$i_h.total,
        type="l",
        lty=1,
        lwd=2,
        col=c("blue","red","darkgreen"),
        xlab="Time (Days)",
        ylab="Number of Infectious",
        main="Baseline vs Intervention")

matlines(time,
         sim_int$epi$i_h.total,
         lty=2,
         lwd=2,
         col=c("blue","red","darkgreen"))

legend("topright",
       legend=c("Baseline mu=0.1",
                "Baseline mu=0.3",
                "Baseline mu=0.5",
                "Intervention mu=0.1",
                "Intervention mu=0.3",
                "Intervention mu=0.5"),
       col=c("blue","red","darkgreen",
             "blue","red","darkgreen"),
       lty=c(1,1,1,2,2,2),
       lwd=2,
       bty="n",
       cex=0.8)

# ======================
# R0 AFTER INTERVENTION (ENHANCED OUTPUT)
# ======================

R0_baseline <- (param$beta_HH * param$mu) / param$rho
R0_intervention <- (beta_HH_int * mu_int) / param$rho

reduction_percent <- 100 * (1 - (R0_intervention / R0_baseline))
delta_R0 <- R0_baseline - R0_intervention

result_table <- data.frame(
  Scenario = c("0.1","0.3","0.5"),
  R0_Baseline = round(R0_baseline,3),
  R0_Intervention = round(R0_intervention,3),
  Reduction_Percent = round(reduction_percent,1),
  Delta_R0 = round(delta_R0,3),
  Baseline_Status = ifelse(R0_baseline > 1, "Supercritical", "Subcritical"),
  Intervention_Status = ifelse(R0_intervention > 1, "Supercritical", "Subcritical")
)

print(result_table)

# ======================
# R0 COMPARISON PLOT
# ======================

R0_baseline <- (param$beta_HH * param$mu) / param$rho
R0_intervention <- (beta_HH_int * mu_int) / param$rho

par(bty="l")
par(mar=c(5,5,4,2))

plot(param$mu, R0_baseline,
     type="b",
     pch=19,
     lwd=2,
     col="blue",
     ylim=c(0, max(R0_baseline)+0.5),
     xlab="mu",
     ylab=expression(R[0]),
     main="Comparison of R0: Baseline vs Intervention")

lines(param$mu, R0_intervention,
      type="b",
      pch=17,
      lwd=2,
      col="darkorange")

abline(h=1, lty=2, lwd=2, col="red")

legend("topright",
       legend=c("Baseline","Intervention","Threshold (R0=1)"),
       col=c("blue","darkorange","red"),
       pch=c(19,17,NA),
       lty=c(1,1,2),
       lwd=2,
       bty="n")

############################################################################

# ======================
# SENSITIVITY HEATMAP
# ======================

# Range parameter
mu_range <- seq(0, 0.6, by = 0.01)
beta_range <- seq(0, 1.2, by = 0.02)

rho_val <- param$rho[1]

# Buat grid kombinasi
R0_matrix <- outer(mu_range, beta_range,
                   function(mu, beta) (beta * mu) / rho_val)

# Plot heatmap
par(mar=c(5,5,4,2))
image(beta_range, mu_range, R0_matrix,
      col = colorRampPalette(c("darkgreen","yellow","red"))(100),
      xlab=expression(beta[HH]),
      ylab=expression(mu),
      main=expression("Sensitivity Heatmap of "*R[0]*" ("*mu*" x "*beta[HH]*")"))

# Tambahkan contour R0 = 1 (threshold)
contour(beta_range, mu_range, R0_matrix,
        levels=1,
        add=TRUE,
        lwd=3,
        col="blue")

# Tandai baseline parameter Anda
points(param$beta_HH[1], param$mu[1], pch=19, col="black")
points(param$beta_HH[1], param$mu[2], pch=19, col="black")
points(param$beta_HH[1], param$mu[3], pch=19, col="black")

legend("topright",
       legend=c("R0 < 1 (Controlled)",
                "R0 > 1 (Outbreak)",
                "R0 = 1 Threshold",
                "Baseline Scenarios"),
       fill=c("darkgreen","red",NA,NA),
       border=NA,
       lwd=c(NA,NA,3,NA),
       col=c(NA,NA,"blue","black"),
       pch=c(NA,NA,NA,19),
       bty="n",
       cex=0.8)

############################################################################