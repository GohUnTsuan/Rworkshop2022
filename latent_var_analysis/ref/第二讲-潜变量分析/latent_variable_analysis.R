#探索性因子分析
library(psych)
library(tidyLPA)#用这个包里面的empathy数据
data.efa<-empathy
data.efa<-data.efa[,-13]
View(data.efa)
dim(data.efa)#467个样本，12个指标

#第一步：利用KMO统计量和Bartlett球形检验判断数据是否适合做因子分析。
correlations <- cor(data.efa) #计算变量相关系数矩阵并赋值给correlations
KMO(correlations) #Overall MSA = 0.86 > 0.7,说明12个指标之间存在较好的相关性
cortest.bartlett(correlations,n=467)# P小于0.05，说明12个指标存在较好的相关性，适合做因子分析
cor.plot(correlations) #可视化相关系数矩阵

#第二步：确定公共因子个数
fa.parallel(correlations, n.obs = 467, fa = "both", n.iter = 100, main = "平行分析碎石图")
#平行分析碎石图建议因子数为3

#第二步：提取因子变量
fa(correlations, nfactors = 3, n.obs = 467, rotate = "none", scores = T, fm = "pa")

#r:相关系数矩阵；
#nfactors:设定提取的因子数，默认为1；
#n.obs:样本量；
#rotate:旋转方法，默认为变异数最小法；
#scores:设定是否计算因子得分，默认不计算；
#fm:设定因子化方法，包含最大似然法(ml)，主轴迭代法(pa)，加权最小二乘法(wls)，广义加权最小二乘法(gls)以及默认的极小残差法(minres)。


#第三步：利用正交/斜交旋转方法使因子变量更具有可解释性
#用正交旋转提取因子
#正交旋转将人为地强制3个因子不相关
fa.varimax <- fa(correlations, nfactors = 3, n.obs = 467, rotate = "varimax", fm = "pa", scores = T)
fa.varimax

#第四步：绘制正交旋转后结果图形
factor.plot(fa.varimax,labels = rownames(fa.varimax$loadings))
fa.diagram(fa.varimax, digits = 3)#这里的系数是因子载荷，不是得分系数

#第五步：计算因子得分
fa.varimax$weights #区分因子载荷和得分系数
#比如PA3 = 0.79*pt5 + 0.20*pt6
#总empathy得分 = a*PA1 + b*PA2 + c*PA3； a,b,c为三个因子对数据的贡献程度百分比。


#验证性因子分析
#安装和加载包
library(semTools)
library(lavaan)
library(semPlot)



#查看lavaan包的示例数据
View(HolzingerSwineford1939)
mydata<-HolzingerSwineford1939

#验证这些观察变量是否能够很好地表示潜变量
model.cfa1<- ' 
         visual  =~ x1 + x2 + x3
	       textual =~ x4 + x5 + x6
	       speed   =~ x7 + x8 + x9 '
fit.cfa1<-lavaan::cfa(model.cfa1,data=mydata)
summary(fit.cfa1,fit.measures=T)
semPaths(fit.cfa1,nCharNodes = 0,whatLabels = "stand",
         edge.color="black",layout = 'tree2',style = "lisrel")


model.cfa2<- ' 
         visual  =~ x1 + x2 + x3
	       textual =~ x4 + x5 + x6 + x3
	       speed   =~ x7 + x8 + x9 '
fit.cfa2<-lavaan::cfa(model.cfa2,data=mydata)
summary(fit.cfa2,fit.measures=T) #基于理论进行了一些调整，貌似好了一点点
semPaths(fit.cfa2,nCharNodes = 0,whatLabels = "stand",
         edge.color="black",layout = 'tree2',style = "lisrel")




#结构方程模型
#下面探索潜变量与潜变量之间的关系
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


#再进阶：潜在增长模型(LGM)——纵向数据的探索
library("semTools")
library("lavaan")
library("semPlot")
library("tidyverse")
library("haven") #import "SPSS data"

#读取数据
data<-haven::read_sav("D:\\somedata\\crime trend data for latent variable analysis\\crime_longitudinal.sav")
#View(data)
mydata<-data%>%select(JANFEB95,MARAPR95,MAYJUN95,JLYAUG95,PA,CPV12590)%>%
  rename(Time1=JANFEB95,Time2=MARAPR95,Time3=MAYJUN95,Time4=JLYAUG95,State=PA,Poverty=CPV12590)%>%
  drop_na()
View(mydata)
dim(mydata)
#可视化先看看前10个被试在四个时期的犯罪倾向走势，对数据有整体的感观
library("ggplot2")
library("ggthemes")
#这里gather()把宽数据转成长数据，方便可视化
data_plot<-mydata%>%mutate(ID=1:952)%>%select(-c("State","Poverty"))%>%
  gather(key=Times,value=Values,-ID)%>%filter(ID<=10)
#可视化
ggplot(data_plot,aes(x=Times,y=Values))+
  geom_line(aes(group=ID),color="black")+ylim(c(2,8))+
  theme_few()
#看看数据整体趋势
data_plot2<-mydata%>%mutate(ID=1:952)%>%select(-c("State","Poverty"))%>%
  gather(key=Times,value=Values,-ID)%>%filter(ID<=50)

data_plot2$Times<-factor(data_plot2$Times)
data_plot2$Times<-as.numeric(data_plot2$Times)
p1<-ggplot(data_plot2,aes(Times,Values))+
  geom_line(aes(group=ID),color="grey",alpha=0.5,size=0.8)+
  geom_smooth(method = lm,formula=y~poly(x,1),color="black",size=1,se=F)+
  geom_smooth(method=lm,formula = y~poly(x,2),color="red",size=1,se=F)+
  geom_smooth(method =lm, formula = y~splines::bs(x, 3),color="orange",size=1,se=F)+
  theme_few()
p1

#下面开始探索潜在增长模型

#首先找到mydata数据集中所有包含“Time”的变量名称
#allTimes<-c("Time1","Time2","Time3","Time4")
allTimes<-grep("Time",names(mydata),value=T)
#设定残差相关
residualTime<-paste("r*",allTimes,sep='')
residual_variance<-paste(allTimes,residualTime,sep="~~",collapse="\n ")
#设定截距负荷
Time_intercept_weight<-paste("1*",allTimes,sep='')
intercept_weight_collapse<-paste(Time_intercept_weight,collapse=" + ")
#设定斜率负荷
Time_slope_weight<-paste(0:(length(allTimes)-1),"*",allTimes,sep='')
slope_weight_collapse<-paste(Time_slope_weight,collapse=" + ")
#设定二次项（用于后面的非线性增长曲线分析）
Time_quadratic_weight<-paste(c(0,1,4,9),"*",allTimes,sep='')
quadratic_weight_collapse<-paste(Time_quadratic_weight,collapse=" + ")

#然后先探索有约束残差的线性潜在增长模型（不包括二次项）
#设立模型
model1<-paste(
  'i=~',intercept_weight_collapse,'\n\n',
  's=~',slope_weight_collapse,'\n\n',
  residual_variance
)
cat(model1)
#分析模型
#用growth函数
fit1<-growth(model1,data=mydata,se="bootstrap",bootstrap=100)
semPaths(fit1)
summary(fit1,fit.measures=F,standardize=F,rsquare=T)
fitmeasures(fit1,c("cfi","tli","gfi","rmsea","srmr"))#fit of goodness
semPaths(fit1,nCharNodes = 0,whatLabels = "standard",
         edge.color="black",layout = 'tree',style = "lisrel")
#####结果发现拟合优度很不错，模型拟合良好，截距因子和斜率因子负相关(-0.023,p<0.01)，
#####说明个体初始犯罪倾向越高，后期犯罪倾向上升越慢
#####截距因子(5.181)和斜率因子(0.113)的方差为0.522和0.017，均有p<0.01，说明个体在最初的犯罪倾向和变化趋势上都有显著差异


#再看看无约束残差的线性潜在增长模型（不包括二次项）
#设立模型
model2<-paste(
  'i=~',intercept_weight_collapse,'\n\n',
  's=~',slope_weight_collapse,'\n\n')
cat(model2)
#分析模型
fit2<-growth(model2,data=mydata,se="bootstrap",bootstrap=100)
semPaths(fit2)
summary(fit2,fit.measures=F,standardize=F,rsquare=T)
fitmeasures(fit2,c("cfi","tli","gfi","rmsea","srmr"))
semPaths(fit2,nCharNodes = 0,whatLabels = "est",
         edge.color="black",layout = 'tree',style = "lisrel")


#再看看无约束条件的非线性增长模型
#设立模型
model3<-paste(
  'i=~',intercept_weight_collapse,'\n\n',
  's=~',slope_weight_collapse,'\n\n',
  'q=~',quadratic_weight_collapse,sep='')
cat(model3)
#分析模型
fit3<-growth(model3,data=mydata,se="bootstrap",bootstrap=100)
semPaths(fit3)
summary(fit3,fit.measures=F,standardize=F,rsquare=T)
fitmeasures(fit3,c("cfi","tli","gfi","rmsea","srmr"))
semPaths(fit3,nCharNodes = 0,whatLabels = "est",
         edge.color="black",layout = 'tree',style = "lisrel")
#模型拟合良好，可以看出个体初始状态与增长速率不相关，二次增长趋势不太明显


#最后看看带协变量的线性增长模型
#比如个体所成长的地区和自身的贫困状态会不会影响犯罪倾向？？？
mydata$interaction<-mydata$State*mydata$Poverty
model4<-paste(
  'i=~',intercept_weight_collapse,'\n\n',
  's=~',slope_weight_collapse,'\n\n',
  'i+s ~ State + Poverty + interaction')
cat(model4)

fit4<-growth(model4,data=mydata,se="bootstrap",bootstrap=100)
semPaths(fit4)
summary(fit4,fit.measures=F,standardize=F,rsquare=T)
fitmeasures(fit4,c("cfi","tli","gfi","rmsea","srmr"))
semPaths(fit4,nCharNodes = 0,whatLabels = "est",
         edge.color="black",layout = 'tree',style = "lisrel")
#####结果发现个体所在的州、本身贫困状态以及二者交互会对个体初始犯罪倾向有显著影响，
#####但是犯罪倾向的变化率(发展趋势)只会收到个体所在州的显著影响。


#那带协变量的非线性增长模型呢？？？自己尝试去解释结果
model5<-paste(
  'i=~',intercept_weight_collapse,'\n\n',
  's=~',slope_weight_collapse,'\n\n',
  'q=~',quadratic_weight_collapse,sep='','\n\n',
  'i+s+q ~ State + Poverty + interaction')
cat(model5)
fit5<-growth(model5,data=mydata,se="bootstrap",bootstrap=100)
semPaths(fit5)
summary(fit5,fit.measures=F,standardize=F,rsquare=T)
fitmeasures(fit5,c("cfi","tli","gfi","rmsea","srmr"))


#潜在剖面模型
library(tidyLPA)
library(tidyverse)
library(cowplot)
library(ggthemes)
#tidyPLA中自带数据集pisaUSA15
View(pisaUSA15)
#dim(pisaUSA15)
data<-pisaUSA15[1:200,]

model <- data %>%
  single_imputation() %>%
  estimate_profiles(n_profiles =1:6)#估计1-6个潜在剖面的结果
model

#或者利用get_fit潜在剖面模型的拟合指数结果
fit_goodness<-get_fit(model)
fit_goodness
#绘制AIC,BIC曲线
p1 <- fit_goodness %>% 
  select(Classes, AIC, BIC) %>% 
  gather("ind","values",-Classes) %>% 
  ggplot(aes(x = Classes, y = values, color = ind, group = ind)) +
  geom_point() +
  geom_line(size = 1) +
  theme_classic()+
  labs(y = "AIC/BIC Values", x = "Classes")
p1
#综合据AIC和BIC的结果，发现class=5/6最好

#绘制熵曲线
#Entropy指的是某种事件不确定性的大小, 值越大表明事件的不确定性越大
p2 <- fit_goodness %>% 
  select(Classes, Entropy) %>% 
  ggplot(aes(x = Classes, y= Entropy)) +
  geom_point() + 
  geom_line(color = "blue",size = 1) +
  theme_classic()
p2
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

#那就分别看看剖面为2/5/6的情况，即2/5/6个剖面(类别)在不同条目上的得分情况
##选择剖面为2的分析结果
#使用get_data函数将上一步得到的类别结果加入到原始数据中
d2 <- get_data(model) %>% 
  filter(classes_number == 2)
names(d2)
#将结果可视化
d2 %>% 
  select(broad_interest, enjoyment, self_efficacy,instrumental_mot, Class) %>%
  mutate(Z1 = (broad_interest-mean(broad_interest))/sd(broad_interest),
         Z2 = (enjoyment-mean(enjoyment))/sd(enjoyment),
         Z3 = (instrumental_mot-mean(instrumental_mot))/sd(instrumental_mot),
         Z4 = (self_efficacy-mean(self_efficacy))/sd(self_efficacy)) %>% 
  group_by(Class) %>%  #将数据分组标准化,去掉测量变量不同量纲的影响
  summarise(mean1 = mean(Z1),
            mean2 = mean(Z2),
            mean3 = mean(Z3),
            mean4 = mean(Z4)) %>% #将标准化的数据放在一列中，以便为后期的分类做准备
  gather("Items","Values",-Class) %>% #把宽数据型转换成长数据型，以便画图
  ggplot(aes(x = as.factor(Class), y = Values, fill = Items, group = Items)) + 
  geom_col(position = "dodge") + #bar chart
  ylab("Z-score") + 
  xlab("Class") + 
  scale_fill_discrete(
    name="Items",
    breaks=c("mean1", "mean2", "mean3","mean4"),
    labels=c("broad_interest", "enjoyment", "self_efficacy","instrumental_mot")) + 
  theme_classic()#修改图例的名称


#换一种表现形式看剖面为5的结果
d5 <- get_data(model) %>% 
  filter(classes_number == 5)
d5<-d5 %>% 
  select(broad_interest, enjoyment, self_efficacy,instrumental_mot, Class) %>%
  mutate(Z1 = (broad_interest-mean(broad_interest))/sd(broad_interest),
         Z2 = (enjoyment-mean(enjoyment))/sd(enjoyment),
         Z3 = (instrumental_mot-mean(instrumental_mot))/sd(instrumental_mot),
         Z4 = (self_efficacy-mean(self_efficacy))/sd(self_efficacy)) %>% 
  group_by(Class) %>% 
  summarise(broad_interest = mean(Z1),
            enjoyment = mean(Z2),
            self_efficacy = mean(Z3),
            instrumental_mot = mean(Z4)) %>% 
  gather("Items", "Zscore",-Class)

d5$Class <- as.factor(d5$Class)

levels(d5$Class) <- c("profile1", "profile2", "profile3", "profile4", "profile5")

ggplot(d5,aes(x = Items, y = Zscore, group = Class)) +
  geom_point(aes(color=Class)) + 
  geom_line(aes(linetype=Class, color=Class),size = 1) +
  scale_linetype_manual(values=c("dotdash","solid", "longdash","dotted","blank"))+
  scale_color_manual(values=c('black','blue','red',"dark green","orange")) +
  theme_classic()

#关于SEM模型被试量N的问题，请查询semPower包