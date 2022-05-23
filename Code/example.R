# p.0
J <- c(0.162, 0.218, 0.142, 0.098, 0.299, 0.561)
# p.1
K <- c(0.283, 0.343, 0.219, 0.083, 0.212, 0.390)
# predicted probabilities
L <- J - K
# observed outcomes
M <- c(0, -1, 1, 0, 1, 0)

plot(L, M)
loess.calibrate <- stats::loess(M ~ L)
N <- predict(loess.calibrate, newdata=L)
print(data.frame(L, M, round(N, 0)))

OC <- abs(mean(M) - mean(N))
ICI <- mean(abs(L - N))
E50 <- median(abs(L - N))
E90 <- as.numeric(quantile(abs(L - N), 0.9))
cat(' OC  :', round(OC, 3), '\n',
    'ICI :', round(ICI, 3), '\n',
    'E50 :', round(E50, 3), '\n',
    'E90 :', round(E90, 3), '\n')

log.loss <- -(as.numeric(M==1)%*%log((1-J)*K)+as.numeric(M==0)%*%log((1-J)*(1-K)+J*K)+as.numeric(M==-1)%*%log(J*(1-K)))/length(L)
Brier <- 1/(2*length(L))*sum(((1-J)*K-as.numeric(M==1))^2+((1-J)*(1-K)+J*K-as.numeric(M==0))^2+(J*(1-K)-as.numeric(M==-1))^2)
cat(' Log loss   :', round(log.loss, 3), '\n',
    'Brier score:', round(Brier, 3), '\n')
