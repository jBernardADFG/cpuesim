# library(readxl)
# dat <- as.data.frame(read_excel("S:/Jordy/louiseOP2020/Data/CPUE.xls", sheet = "CPUEL", col_names = F))
# colnames(dat) <- c("tran", "set", "ct")
# dat[,1][dat[,1]==2] <- 20
# dat <- dat[order(dat$tran, dat$set),]
# tran <- unique(dat[,1])
# df <- data.frame()
# for(i_1 in 1:nrow(dat)){
#   for(i_2 in 1:nrow(dat)){
#     if (dat[i_1,1]==dat[i_2,1]){ # sets are in the same transect
#       d <- abs(dat[i_1,2]-dat[i_2,2])
#       d <- d*125
#       v <- (dat[i_1,3]-dat[i_2,3])^2
#       df <- rbind(df, c(d,v))
#     }
#   }
# }
# df <- df[order(df[,1]),]
# df <- df[!df[,1]==0,]
# rdf <- data.frame()
# for(d in unique(df[,1])){
#   ddf <- df[df[,1]==d,]
#   sv <- sum(ddf[,2])/(nrow(ddf)) # Empirical semivariance
#   n <- nrow(ddf)
#   rdf <- rbind(rdf, c(d, sv, n))
# }
# rdf <- data.frame(rdf)
# rdf <- rdf[rdf[,1]<3000,]
# names(rdf) <- c("gam", "dist", "np")
# gam <- rdf[,2]
# dist <- rdf[,1]
# weights <- rdf[,3]
# nug <- 1.2
# y <- gam-nug
# r_grid <- 1:3000
# coeffs <- logliks <- rep(NA, length(r_grid))
# vf <- varFunc(~weights)
# for (i in 1:length(r_grid)){
#   r <- r_grid[i]
#   x <- 1-exp(-dist/r)
#   fm1 <- gls(y ~ 0 + x, weights = vf)
#   coeffs[i] <- fm1$coefficients[1]
#   logliks[i] <- fm1$logLik
# }
# best_fit <- which(logliks == max(logliks))
# p_sill <- coeffs[best_fit]
# r <- r_grid[best_fit]
# x_grid <- 1:3000
# y_grid <- p_sill*(1-exp(-x_grid/r))+nug
# plot(0, bty = 'n', pch = '', ylab = "semivariance", xlab = "distance (m)", xlim=c(0, 3000), ylim=c(1.2, 1.9))
# with(rdf, symbols(x=rdf[,1], y=rdf[,2], circles=rdf[,3], inches=1/5,
#                   ann=F, bg="steelblue2", fg=NULL, add=T)) # Maybe add some labels "dist (m)" and "semivariance"
# lines(x_grid, y_grid)
# # Let's see if we can get a correlogram 
# y_grid <- (nug+p_sill-y_grid)/(nug+p_sill)
# plot(x_grid, y_grid, ty="l", x_lab="distance (m)", y_lab="correlation") # Here's our correlogram
# dat <- data.frame(distance=x_grid, correlation=y_grid)
# # Basic line plot
# library(ggplot2)
# ggplot(data=dat, aes(x=distance, y=correlation))+
#   geom_line()+ ggtitle("Correlogram") +
#   xlab("distance (m)") + ylab("correlation")+theme(plot.title = element_text(hjust=0.5, size=14, face="bold.italic"),
#                      axis.title.x = element_text(size=14, face="bold"),
#                      axis.title.y = element_text(size=14, face="bold"))
