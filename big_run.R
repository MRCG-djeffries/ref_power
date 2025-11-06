this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(file.path(this_dir, "run_sim3.R")) 
big_run = function(){
nrow = 10
vac = pow = matrix(data =0, nrow =nrow, ncol = 11)
Nevents = 20:30
for (i in 1 : 11){
  for ( j in 1 : nrow){
      x1=run_sim3(nevents = Nevents[i],looks_at = (i-1) + c(10, 15, 20),rho_target = NULL)
      pow[j,i]=x1$power
      vac[j,i]=x1$meanVacc
  }
}
saveRDS(pow,"powmat3.rds")
saveRDS(vac,"vacmat3.rds")

nrow = 10
vac = pow = matrix(data =0, nrow =nrow, ncol = 11)
Nevents = 20:30
for (i in 1 : 11){
  for ( j in 1 : nrow){
    x1=run_sim3(nevents = Nevents[i],looks_at = (i-1) + c(10, 15, 20),rho_target = 0.3)
    pow[j,i]=x1$power
    vac[j,i]=x1$meanVacc
  }
}
saveRDS(pow,"powmat4.rds")
saveRDS(vac,"vacmat4.rds")

# Read simulation results
pow1 <- readRDS("powmat3.rds")   # no frailty
vac1 <- readRDS("vacmat3.rds")
pow2 <- readRDS("powmat4.rds")   # latent frailty
vac2 <- readRDS("vacmat4.rds")

x  <- 20:30
p1 <- 100 * colMeans(pow1)
p2 <- 100 * colMeans(pow2)
v1 <- colMeans(vac1)
v2 <- colMeans(vac2)
par(mar=c(4,4,1,1))
# --- Plot 1: Power ---
plot(loess.smooth(x, p1), type = "l", lwd = 3, col = "blue",
     xlab = "Target number of cases",
     ylab = "Power (%)",
     ylim = c(75, 100))
#main = "Simulated Power by Event Target")
lines(loess.smooth(x, p2), col = "red", lwd = 3)
grid()
legend("topleft",
       legend = c("No frailty", "Latent frailty"),
       col = c("blue", "red"), lwd = 3, bty = "n")

# --- Plot 2: Vaccinated subjects enrolled ---
plot(loess.smooth(x, v1), type = "l", lwd = 3, col = "blue",
     xlab = "Target number of cases",
     ylab = "Mean number vaccinated at endpoint",
     ylim = c(100, 160))
#main = "Vaccinated Subjects by Event Target")
lines(loess.smooth(x, v2), col = "red", lwd = 3)
grid()
legend("topleft",
       legend = c("No frailty", "Latent frailty"),
       col = c("blue", "red"), lwd = 3, bty = "n")


}