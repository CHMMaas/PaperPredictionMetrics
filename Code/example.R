n <- 8

# matched.p.0
J <- c(0.162, 0.218, 0.142, 0.098, 0.299, 0.561, 0.243, 0.345)
# J <- runif(n)
# matched.p.1
K <- c(0.283, 0.343, 0.219, 0.083, 0.212, 0.390, 0.201, 0.199)
# predicted probabilities, matched.tau.hat
L <- J - K
# observed outcomes, matched.tau.obs
M <- c(0, -1, 1, 0, 1, 0, 0, 1)
# M <- sample(-1:1, n, replace=TRUE)

plot(L, M)
loess.calibrate <- stats::loess(M ~ L)
N <- predict(loess.calibrate, newdata=L, type="response")

OC <- abs(mean(M) - mean(N))
ICI <- mean(abs(L - N))
E50 <- median(abs(L - N))
E90 <- as.numeric(quantile(abs(L - N), 0.9))
n.p <- length(L)
t.1 <- (1-K)*J
t.0 <- (1-K)*(1-J)+K*J
t.min1 <- K*(1-J)
I.1 <- M==1
I.0 <- M==0
I.min1 <- M==-1
log.loss <- -(as.numeric(M==1)%*%log((1-K)*J)+
                as.numeric(M==0)%*%log((1-K)*(1-J)+K*J)+
                as.numeric(M==-1)%*%log(K*(1-J)))/n.p
Brier <- 1/(2*n.p)*sum(((1-K)*J-as.numeric(M==1))^2+
                         ((1-K)*(1-J)+K*J-as.numeric(M==0))^2+
                         (K*(1-J)-as.numeric(M==-1))^2)

matched.patients <- data.frame(subclass=1:n, matched.tau.obs=M, 
                               matched.p.0=J, matched.p.1=K, matched.tau.hat=L)
EB <- HTEPredictionMetrics::E.for.Benefit(matched.patients=matched.patients,
                                    CI=FALSE, message=FALSE, replace=FALSE)
OP <- HTEPredictionMetrics::OP.for.Benefit(matched.patients=matched.patients,
                                     CI=FALSE, message=FALSE, replace=FALSE)

# manual
cat(' OC  :', round(OC, 3), '\n',
    'ICI :', round(ICI, 3), '\n',
    'E50 :', round(E50, 3), '\n',
    'E90 :', round(E90, 3), '\n',
    'Log loss   :', round(log.loss, 3), '\n',
    'Brier score:', round(Brier, 3), '\n')

# package
cat(' OC  :', round(OC, 3), '\n',
    'ICI :', round(EB$Eavg.for.benefit, 3), '\n',
    'E50 :', round(EB$E50.for.benefit, 3), '\n',
    'E90 :', round(EB$E90.for.benefit, 3), '\n',
    'Cross entropy:', round(OP$Cross.entropy.for.benefit, 3), '\n',
    'Brier score  :', round(OP$Brier.for.benefit, 3), '\n')

print(round(data.frame(J, K, L, M, N), 3))
