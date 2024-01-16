Cadenas_sim_moderada <- list("B_POST" = B_POST, 
                             "LV_POST" = LV_POST,
                             "NU_POST" = NU_POST, 
                             "PI_POST" = PI_POST,
                             "SIGMA_POST" = SIGMA_POST,
                             "TAU_POST" = TAU_POST, 
                             "PHI_REAL" = Phi, 
                             "PI_REAL" = rep(c(PI_1),3))

Cadenas_sim_extrema <- list("B_POST" = B_POST, 
                            "LV_POST" = LV_POST,
                            "NU_POST" = NU_POST, 
                            "PI_POST" = PI_POST,
                            "SIGMA_POST" = SIGMA_POST,
                            "TAU_POST" = TAU_POST, 
                            "PHI_REAL" = Phi, 
                            "PI_REAL" = rep(c(PI),3))


Cadenas_sim_corta   <- list("B_POST" = B_POST, 
                            "LV_POST" = LV_POST,
                            "NU_POST" = NU_POST, 
                            "PI_POST" = PI_POST,
                            "SIGMA_POST" = SIGMA_POST,
                            "TAU_POST" = TAU_POST, 
                            "PHI_REAL" = Phi, 
                            "PI_REAL" = rep(c(PI),3))

save(list=c("Cadenas_sim_moderada", 
            "Cadenas_sim_extrema",
            "Cadenas_sim_corta"), file = "Cadenas_Consolidadas.Rdata")


.B_POST  = Cadenas_sim_extrema[["B_POST"]] |> Reduce(f = 'rbind')
.PI_POST = Cadenas_sim_extrema[["PI_POST"]] |> Reduce(f = 'rbind')

comp_phi <- .B_POST[,!(seq_along(b_list[[1]]) %% (k[1]+r) %in% c(1:r))]

par(mfrow = c(6,8), mar=c(1,1,1,0))
for (i in 1:48) {
  comp_phi[,i] |> hist(breaks = 50, main = NULL, xlim = c(-1,1))
  abline(v = Cadenas_sim_corta[['PHI_REAL']][i], col='red', lwd=3)
  #abline(v = mean(comp_phi[,i]), col = 'blue', lwd = 3)
  abline(v = quantile(comp_phi[,i], c(0.025, 0.975)), col = 'green', lwd = 3)
}

comp_pi <- .PI_POST
par(mfrow = c(6,8), mar=c(1,1,1,0))
for (i in 1:48) {
  comp_pi[,i] |> hist(breaks = 50, main = NULL, xlim = quantile(comp_pi[,i], c(0.025, 0.975)))
  abline(v = Cadenas_sim_corta[['PI_REAL']][i], col='red', lwd=3)
  #abline(v = mean(comp_pi[,i]), col = 'blue', lwd = 3)
  abline(v = quantile(comp_pi[,i], c(0.025, 0.975)), col = 'green', lwd = 3)
}
