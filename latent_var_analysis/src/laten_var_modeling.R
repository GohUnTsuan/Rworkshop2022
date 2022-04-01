################################################################################
############### Tsinghua Computional Social Science R Workshop 2022 Spring
############### Latent Variable Modeling
############### Wenquan Wu
############### 2022-04-01 
################################################################################
# load ----
library(tidyverse)
load('latent_var_analysis/data/wvs_wave7.rdata')
data <- `WVS_Cross-National_Wave_7_R_v3_0`

dat_eco <- data %>% 
  select(
    Q33,
    Q34,
    Q35,
    Q36,
    Q37,
    Q38,
    Q39,
    Q40,
    Q41
  )

dat_crime <- data %>% 
  select(
    Q177,
    Q178,
    Q179,
    Q180,
    Q181,
    Q182,
    Q183,
    Q184,
    Q185,
    Q186,
    Q187,
    Q188,
    Q189,
    Q190,
    Q191,
    Q192,
    Q193,
    Q194,
    Q195
  )
dat_crime[dat_crime == -5] <- NA_real_
dat_crime[dat_crime == -4] <- NA_real_
dat_crime[dat_crime == -3] <- NA_real_
dat_crime[dat_crime == -2] <- NA_real_
dat_crime[dat_crime == -1] <- NA_real_

dat_eco <- dat_crime
table(dat_eco$Q33)
table(dat_eco$Q34)
table(dat_eco$Q36)

dat_eco[dat_eco == -5] <- NA_real_
dat_eco[dat_eco == -2] <- NA_real_
dat_eco[dat_eco == -1] <- NA_real_
dat_eco[dat_eco == -4] <- NA_real_
# Principal components analysis (PCA) ----
#探索性因子分析
library(psych)
library(tidyLPA)#用这个包里面的empathy数据
library(pacman)
p_load(
  tidyLPA,
  corrplot
)
data.efa<-dat_eco %>% na.omit()
# data.efa<-empathy
# data.efa<-data.efa[,-13]
dim(data.efa) #78035个样本，12个指标
n <- dim(data.efa)[1]
res <- cor.mtest(data.efa, conf.level = 0.95)

#第一步：利用KMO统计量和Bartlett球形检验判断数据是否适合做因子分析。
correlations <- cor(data.efa) #计算变量相关系数矩阵并赋值给correlations
KMO(correlations) # Overall MSA = 0.86 > 0.7,说明12个指标之间存在较好的相关性
print(KMO(correlations)$Image) # Image
cortest.bartlett(correlations,n=n)# P小于0.05，说明12个指标存在较好的相关性，适合做因子分析
# cor.plot(correlations) #可视化相关系数矩阵
corrplot(correlations, method = "ellipse", type = "lower", p.mat = res$p, sig.level = 0.05, order = "hclust")

cortest.bartlett(correlations, n=n) # bartlett test
alpha(data.efa) # cronbach alpha

#第二步：确定公共因子个数
fa.parallel(correlations, n.obs = n, fa = "both", n.iter = 100, main = "平行分析碎石图")
#平行分析碎石图建议因子数为3

#第二步：提取因子变量
fa(correlations, nfactors = 3, n.obs = n, rotate = "none", scores = T, fm = "pa")

#r:相关系数矩阵；
#nfactors:设定提取的因子数，默认为1；
#n.obs:样本量；
#rotate:旋转方法，默认为变异数最小法；
#scores:设定是否计算因子得分，默认不计算；
#fm:设定因子化方法，包含最大似然法(ml)，主轴迭代法(pa)，加权最小二乘法(wls)，广义加权最小二乘法(gls)以及默认的极小残差法(minres)。


#第三步：利用正交/斜交旋转方法使因子变量更具有可解释性
#用正交旋转提取因子
#正交旋转将人为地强制3个因子不相关
fa.varimax <- fa(correlations, nfactors = 2, n.obs = n, rotate = "varimax", fm = "pa", scores = T)
fa.varimax

#第四步：绘制正交旋转后结果图形
factor.plot(fa.varimax,labels = rownames(fa.varimax$loadings))
fa.diagram(fa.varimax, digits = 3)#这里的系数是因子载荷，不是得分系数

#第五步：计算因子得分
fa.varimax$weights #区分因子载荷和得分系数
#比如PA3 = 0.79*pt5 + 0.20*pt6
#总empathy得分 = a*PA1 + b*PA2 + c*PA3； a,b,c为三个因子对数据的贡献程度百分比  

# Item Response Theory (IRT) ----
# Load the mirt library
p_load(mirt)

# Fit
data.irt <- data.efa
mod1 <- (mirt(data.irt, 2, verbose = FALSE, itemtype = 'graded', SE = TRUE))
M2(mod1, type = "C2", calcNULL = FALSE) 
# RMSEA < 0.06 SRMSR < 0.08 TLI > 0.9 CFI > 0.9

# Item fit
itemfit(mod1)
RMSEA.S_X2 < 0.06

# IRT para
IRT_parms <- coef(mod1, IRTpars = TRUE, simplify = TRUE)
IRT_parms$items

# Factor loadings
summary(mod1)
# loading > 0.5

# Category characteristic curves
plot(mod1, type='trace', which.item = c(1,2,3,4,5,6), facet_items=T, 
     as.table = TRUE, auto.key=list(points=F, lines=T, columns=4, space = 'top', cex = .8), 
     theta_lim = c(-3, 3), 
     main = "")
# 效果爆炸，重新做undimension

# factor 1
data.irt1 <- data.irt %>% select(
  Q37,
  Q38,
  Q39,
  Q40,
  Q41
)
mod1 <- (mirt(data.irt1, 1, verbose = FALSE, itemtype = 'graded', SE = TRUE))
M2(mod1, type = "C2", calcNULL = FALSE) 
# RMSEA < 0.06 SRMSR < 0.08 TLI > 0.9 CFI > 0.9

# Item fit
itemfit(mod1)
# RMSEA.S_X2 < 0.06

# IRT para
IRT_parms <- coef(mod1, IRTpars = TRUE, simplify = TRUE)
IRT_parms$items

# Factor loadings
summary(mod1)
# loading > 0.5

# Category characteristic curves
plot(mod1, type='trace', which.item = c(1,2,3,4), facet_items=T, 
     as.table = TRUE, auto.key=list(points=F, lines=T, columns=4, space = 'top', cex = .8), 
     theta_lim = c(-3, 3), 
     main = "")

# Item information curves
plot(mod1, type='infotrace', which.item = c(1,2,3,4), facet_items=T, 
     as.table = TRUE, auto.key=list(points=F, lines=T, columns=1, space = 'right', cex = .8), 
     theta_lim = c(-3, 3), 
     main="")

# Scale information and conditional standard errors
plot(mod1, type = 'infoSE', theta_lim = c(-3, 3), 
     main="")

# Conditional reliability
plot(mod1, type = 'rxx', theta_lim = c(-3, 3), 
     main="" )
marginal_rxx(mod1)

# Scale characteristic curve
plot(mod1, type = 'score', theta_lim = c(-3, 3), main = "")

# factor 2
data.irt2 <- data.irt %>% select(
  Q33,
  Q34,
  Q35,
  Q36)
mod2 <- (mirt(data.irt2, 1, verbose = FALSE, itemtype = 'graded', SE = TRUE))
M2(mod2, type = "C2", calcNULL = FALSE) 
# RMSEA < 0.06 SRMSR < 0.08 TLI > 0.9 CFI > 0.9

# Item fit
itemfit(mod2)
# RMSEA.S_X2 < 0.06

# IRT para
IRT_parms <- coef(mod2, IRTpars = TRUE, simplify = TRUE)
IRT_parms$items

# Factor loadings
summary(mod2)
# loading > 0.5

# Category characteristic curves
plot(mod2, type='trace', which.item = c(1,2,3,4), facet_items=T, 
     as.table = TRUE, auto.key=list(points=F, lines=T, columns=4, space = 'top', cex = .8), 
     theta_lim = c(-3, 3), 
     main = "")

# Item information curves
plot(mod2, type='infotrace', which.item = c(1,2,3,4), facet_items=T, 
     as.table = TRUE, auto.key=list(points=F, lines=T, columns=1, space = 'right', cex = .8), 
     theta_lim = c(-3, 3), 
     main="")

# Scale information and conditional standard errors
plot(mod2, type = 'infoSE', theta_lim = c(-3, 3), 
     main="")

# Conditional reliability
plot(mod2, type = 'rxx', theta_lim = c(-3, 3), 
     main="" )
marginal_rxx(mod1)

# Scale characteristic curve
plot(mod2, type = 'score', theta_lim = c(-3, 3), main = "")


# Structural Equation Modeling (SEM) ----
#结构方程模型
#下面探索潜变量与潜变量之间的关系
mydata<-HolzingerSwineford1939
model.sem<-'
         visual  =~ x1 + x2 + x3          #外生潜变量(自变量)
	       textual =~ x4 + x5 + x6         #外生潜变量(自变量)
	       speed   =~ x7 + x8 + x9          #内生潜变量(因变量)
         speed   =~ visual + textual'     #结构方程模型
fit.sem<-lavaan::sem(model.sem,data=mydata)
summary(fit.sem,fit.measures=T)
semPaths(fit.sem,nCharNodes = 0,whatLabels = "stand",
         edge.color="black",layout = 'tree2',style = "lisrel")

#进阶：看看两个潜变量的交互作用
#由于A*B源于A和B，所以相关性很高，会存在共线性问题
#利用semTools包对自变量进行中心化处理，从而避免该问题
#进行基于限制中心化策略的公式是：Y=a0'+a1'(X1-X1@)+a2'(X2-X2@)+a3'(X1-X1@)(X2-X2@)+e'
#通常有中心化策略、正交化策略、双中心化策略等，这里推荐使用双中心化策略
data.dmc<-indProd(data=mydata,
                  var1=c("x1","x2","x3"),
                  var2=c("x4","x5","x6"),
                  doubleMC=T)
View(data.dmc)
model.dmc<-'
      visual  =~ x1 + x2 + x3                             #外生潜变量(自变量)
	    textual =~ x4 + x5 + x6                             #外生潜变量(自变量)
	    visual.textual =~ x1.x4 + x2.x5 + x3.x6             #经过双中心化处理的交互潜变量
	    speed   =~ x7 + x8 + x9                             #内生潜变量(因变量)
      speed   =~ visual + textual + visual.textual'       #结构方程模型
#切记这里的交互要写成A.B的形式
fit.dmc<-lavaan::sem(model.dmc,data = data.dmc)
summary(fit.dmc,fit.measures=T)

#可视化结构方程图
semPaths(fit.dmc,nCharNodes = 0,whatLabels = "est",
         edge.color="black",layout = 'tree2',style = "lisrel")

#基于和正交化策略处理共线性问题
data.orth<-orthogonalize(data=mydata,
                         var1=c("x1","x2","x3"),
                         var2=c("x4","x5","x6"))
View(data.orth)
model.orth<-'
      visual  =~ x1 + x2 + x3                             #外生潜变量(自变量)
	    textual =~ x4 + x5 + x6                             #外生潜变量(自变量)
	    visual.textual =~ x1.x4 + x2.x5 + x3.x6             #经过正交化处理的交互潜变量
	    speed   =~ x7 + x8 + x9                             #内生潜变量(因变量)
      speed   =~ visual + textual + visual.textual' 
fit.orth<-lavaan::sem(model.orth,data = data.orth)
summary(fit.orth,fit.measures=T)
semPaths(fit.orth,nCharNodes = 0,whatLabels = "est",
         edge.color="black",layout = 'tree2',style = "lisrel")
#如果比较两个模型
#可以用semTools::compareFit(model1,model2)
comparefit<-semTools::compareFit(fit.dmc,fit.orth)
summary(comparefit)