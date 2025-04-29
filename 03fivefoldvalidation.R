#long-short-survival project
#MacOS 15.4.1
#R 4.4.2
#Jiacheng Dai


library(caret)
library(pROC)

datMeta1 = read.csv("02result/03model-construct/ImmuneModel.longshort.csv",header=T,row.names=1)

vars = c("PFS.norm","PFS.end","main_TNM","PDL1_mRNA")
data = datMeta1[,vars]

set.seed(666)
seed = createFolds(y = data$PFS.end, k=5)

data_test = data[seed[[4]],]
data_train = data[-seed[[4]],]

seed_pre = coxph(Surv(PFS.norm, PFS.end)~main_TNM + PDL1_mRNA, data=data_train)
seed_predict = predict(seed_pre, type = "survival", newdata = data_test)
roc1 = roc((data_test[,2]),seed_predict)
round(auc(roc1),3)
round(ci(roc1),3)
plot(roc1, print.auc=T, auc.polygon=T, grid=c(0.1, 0.2),
		grid.col=c('green','red'), max.auc.polygon=T,
		auc.polygon.col='skyblue',
		print.thres=T)

auc.value = as.numeric()
for (i in 1:5) {
	data_test = data[seed[[i]],]
	data_train = data[-seed[[i]],]
	seed_pre = coxph(Surv(PFS.norm, PFS.end)~main_TNM + PDL1_mRNA, data=data_train)
	seed_predict = predict(seed_pre, type = "survival", newdata = data_test)
	auc.value = append(auc.value, as.numeric(auc(as.numeric(data_test[,2]),seed_predict)))
}
auc.value
#[1] 0.5666667 0.9500000 0.5490196 0.7500000 0.9166667
mean(auc.value)
#[1] 0.7464706



attach(data)
pfs1 = ts(PFS.norm, start=1, frequency=1)
pfs2 = ts(PFS.end, start=1, frequency=1)
tnm1 = ts(main_TNM, start=1, frequency=1)
pdl1 = ts(PDL1_mRNA, start=1, frequency=1)
detach(data)

window_size = 80
step = 1

cv_r <- list()
RMSE <- list()
MAE <- list()

set.seed(8888)

for (i in seq(window_size, length(pfs1), by = step)) {
  train_X <- window(pfs1, end = i - 30)
  train_Y <- window(pfs2, end = i - 30)
  train_Z <- window(tnm1, end = i - 30)
  train_W <- window(pdl1, end = i - 30)
  
  test_X <- window(pfs1, start = i-29, end = i)
  test_Y <- window(pfs2, start = i-29, end = i)
  test_Z <- window(tnm1, start = i-29, end = i)
  test_W <- window(pdl1, start = i-29, end = i)
  
  model <- coxph(Surv(train_X, train_Y) ~ train_Z + train_W)
  
  prediction <- predict(model, newdata = data.frame(train_X = test_X, train_Y = test_Y, train_Z = test_Z, train_W = test_W))
  
  # 选取一些指标来评价模型
  cv_r[[i]] <- summary(model)$coefficient[1,5]
  RMSE[[i]] <- sqrt(mean((na.omit(prediction-test_Y))^ 2))
  MAE[[i]] <- mean(abs(na.omit(prediction-test_Y)))
}

cv_r_values <- unlist(cv_r)
cv_RMSE_values <- unlist(RMSE)
cv_MAE_values <- unlist(MAE)

