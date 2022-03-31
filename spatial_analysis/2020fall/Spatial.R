

#
install.packages("readxl")
#maps
#install.packages("sp")
install.packages("tmap")
#install.packages("leaflet")
#install.packages("mapview")
#回归
#
install.packages("spdep")
install.packages("splm")
#install.packages("rgdal")
#install.packages("spatialreg")

#开始咯
library(readxl)
library("spdep")
head(pdata)
#面板 类似于Stata :xtset i t
library(plm)
datap <- pdata.frame(pdata,c("city","year"))
names(datap)


######我是羞羞的顶刊分界线#######

#时间滞后
lag(datap$gdp,1)
#导入并转换权重矩阵
w1 <- read_xlsx("./spatial_analysis/2020fall/data/w1.xlsx",col_names = FALSE)
w2 <- as.matrix(w1)
w <- mat2listw(w2,style="W")
#OLS
pool = plm(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu, data=datap, model="pooling")
summary(pool)
#个体效应和时间效应
#系数
plmtest(pool, effect="twoways", type="honda")
#FE
FEeffect = plm(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu, data=datap, model="within", effect="individual")
summary(FEeffect)
#or
LSDVindiv = lm(gdp ~ 0 + kj + l + ks + pe + inex + new_inc + pri_en + high_stu + factor(city), data=datap)
summary(LSDVindiv)
#LM test:空间误差or空间滞后
library(splm)
slmtest(FEeffect, listw=w, test="lml")
slmtest(FEeffect, listw=w, test="rlml")
slmtest(FEeffect, listw=w, test="lme")
slmtest(FEeffect, listw=w, test="rlme")

##大头来了 面板空间计量模型估计
#排序：first time, second region
library(splm)
datasp=datap[order(datap$year),]
datasp
#一、 空间权重矩阵

# 导入并合并数据
library(readxl)
library(rgdal)
library(spdep)
 # 导入sha文件
library(rgdal)
# setwd('C:\\Users\\grx\\Desktop\\club_R_Spatial\\data')
shpt <- readOGR("./spatial_analysis/2020fall/data/广东地级市.shp")

 #可视化
plot(shpt)
# 合并经济变量数据
cdatashpt <- merge(shpt, cdata, by = "city")
# 1.spdep 用于生成邻接矩阵
library("spdep")
queen.w <- poly2nb(cdatashpt, row.names = cdatashpt$city, queen = TRUE)
summary(queen.w)
#转换格式 nb → listw
queen.wl <- nb2listw(queen.w, style = "W")
summary(queen.wl)
#rook矩阵
rook.w <- poly2nb(cdatashpt, row.names = cdatashpt$city, queen = FALSE)
summary(rook.w)
#绘制权重矩阵 <仅支持 readOGR 等读取的格式>
queen.w <- poly2nb(cdatashpt, row.names = cdatashpt$city, queen = TRUE)
rook.w <- poly2nb(cdatashpt, row.names = cdatashpt$city, queen = FALSE)
plot(cdatashpt, border = "black")
plot(queen.w, coordinates(cdatashpt), add = TRUE, col = "red", lwd = 2)
plot(rook.w, coordinates(cdatashpt), add = TRUE, col = "yellow")
# 2.k近邻
 # 提取坐标值
coords <- coordinates(cdatashpt)
head(coords, 5)
# if longlat = TRUE, circle distances are used. objects k1neigh and k2neigh are of class knn.
k1neigh <- knearneigh(coords, k = 1, longlat = TRUE) # 1-nearest neighbor
k2neigh <- knearneigh(coords, k = 2, longlat = TRUE) # 2-nearest neighbor
# 3.反距离矩阵
# Inverse weight matrix
dist.mat <- as.matrix(dist(coords, method = "euclidean"))
dist.mat[1:5, 1:5]
# 计算反距离矩阵
dist.mat.inv <- 1 / dist.mat # 1 / d_{ij}
diag(dist.mat.inv) <- 0 # 0 in the diagonal
dist.mat.inv[1:5, 1:5]
# 4.标准化权重矩阵
# Standardized inverse weight matrix
dist.mat.inve <- mat2listw(dist.mat.inv, style = "W", row.names = cdatashpt$city)
summary(dist.mat.inve)
# 5.导入自用的 excel 矩阵
library(readxl)
w1 <- read_xlsx("./spatial_analysis/2020fall/data/w1.xlsx", col_names = FALSE)
# 转化为listw格式
w2 <- as.matrix(w1)
w <- mat2listw(w2, style="W")
w2[1:5, 1:5]

#二、空间自相关检验
#Moran'I 指数
 moran(cdatashpt$gdp2017, listw=w, n=length(cdatashpt$gdp2017), S0=Szero(w))
 #出现孤岛或缺失值
  ##zero.policy默认值为空，使用全局选项值；如果为真，则将0分配给没有邻居的区域的滞后值，如果为假，则分配为NA
  ##NAOK如果 TRUE，则 x 中的任何 NA 或 NaN 或 Inf 值都将传递给外部函数。如果 FALSE ，则存在 NA 或 NaN 或 Inf 值将被视为错误。
 # Monte-Carlo simulation of Moran's I
 set.seed(12345)
 moran.mc(cdatashpt$gdp2017, listw = w, nsim = 999, alternative = 'greater')
 # Moran 散点图
 moran.plot(cdatashpt$gdp2017, w, zero.policy=NULL, spChk=NULL, labels=TRUE, xlab=NULL, ylab=NULL, quiet=NULL)
 # (1) Moran’s I test for spatial autocorrelation
 moran.test(cdatashpt$gdp2017, w)
 # (2) Geary test
 geary.test(cdatashpt$gdp2017, listw=w, randomisation=TRUE, alternative="greater")
 #2. Geary’s C 检验
  # 计算 Geary’s C 统计量
 geary(cdatashpt$gdp2017, listw=w, n=length(w), n1=length(w)-1, S0=Szero(w))
  # Monte-Carlo simulation of Geary's C
  geary.mc(cdatashpt$gdp2017, listw=w, nsim=999, alternative="greater")
  # 模拟分布图
 set.seed(12345)
 gdpgeary <- geary.mc(cdatashpt$gdp2017, listw=w, nsim=999, alternative="greater")
 plot(gdpgeary, type='l', col='orange')
 gdpgeary.dens = density(gdpgeary$res)
 polygon(gdpgeary.dens, col="gray")
 abline(v=gdpgeary$statistic, col='orange',lwd=2)
 #3. Getis-Ord global G 检验
 globalG.test(cdatashpt$gdp2017, listw=w, alternative="greater")
  #空间相关图  
  w.nb <- w$neighbours
  spcorrI = sp.correlogram(w.nb, cdatashpt$gdp2017, order = 2, method = "I", style = "W", randomisation = TRUE)
  spcorrI
  plot(spcorrI, main="Spatial correlogram of gdp2017")
# 局部空间自相关检验
  # 局部moran检验
  localmoran(cdatashpt$gdp2017, listw=w, alternative = "greater")
  # Local Getis-Ord Gi 和 Gi*统计量
  localG(cdatashpt$gdp2017, listw=w)
# 基于残差项的moran检验
 ols <- lm(gdp2017 ~ kj2017 + l2017 + ks2017 + pe2017 + inex2017 + new_inc2017 + pri_en2017 + high_stu2017, data = cdatashpt)
 lm.morantest(ols, listw = w, alternative = "two.sided")
 # 散点图的绘制
 moran.plot(ols$residuals, w)
 # 局部moran检验
 localmoran(ols$residuals, w)
 
# 三、地图Maps
 library(readxl)
 # setwd('C:\\Users\\grx\\Desktop\\club_R_Spatial\\data')
 cdata <- read_xlsx("./spatial_analysis/2020fall/data/cdata.xlsx")
 # 导入地图文件
 library(sf)
 shp <- st_read("./spatial_analysis/2020fall/data/广东地级市.shp")
 plot(shp, main = "广东省地级市", axes = TRUE)
 # 合并数据
 cdatashp <- merge(shp, cdata, by = "city")
 # 空间可视化
 library(tmap)
 library(leaflet)
 library(mapview)
 # 行政图
 ma1 = tm_shape(cdatashp) + tm_fill(col = "city")
 ma1
 # 行政图 + 边界
 ma2 = tm_shape(cdatashp) + tm_fill(col = "city") + tm_borders()
 ma2
 # 分级显示
  ma4 = tm_shape(cdatashp) + tm_polygons(col = "gdp2017")
  ma4
 # 自定义分级显示
 ma5 = tm_shape(cdatashp) + tm_polygons(col = "gdp2017", breaks = c(0, 2, 4, 6, 8, 10))
 ma5
 # 更换主题
 ma4 + tm_style("classic")
 ma4 + tm_style("cobalt")
 ma4 + tm_style("col_blind")
 # 图片布局
 # geom_sf() or geom_
 tmap_arrange(ma1, ma5, ma2, ma4)
 # 添加图例
 ma4 +
   tm_compass(type = "8star", position = c("left", "top")) +
   tm_scale_bar(breaks = c(0, 100, 200), text.size = 1)
 
 library(ggplot2)
 g1 = ggplot() + geom_sf(data = cdatashp, aes(fill = gdp2017)) +
   geom_sf() +
   scale_x_continuous(breaks = c(0, 4, 8))
 g1
 
 # 四、LM test:空间误差or空间滞后
 library(splm)
 slmtest(FEeffect, listw=w, test="lml")
 slmtest(FEeffect, listw=w, test="rlml")
 slmtest(FEeffect, listw=w, test="lme")
 slmtest(FEeffect, listw=w, test="rlme")
 
 # 五、回归分析
 # 空间面板滞后项
datasp$wtx_gdp <- slag(datasp$gdp, w, 1)
datasp$wtx_kj <- slag(datasp$kj, w, 1)
datasp$wtx_l <- slag(datasp$l, w, 1)
datasp$wtx_ks <- slag(datasp$ks, w, 1)
datasp$wtx_pe <- slag(datasp$pe, w, 1)
datasp$wtx_inex <- slag(datasp$inex, w, 1)
datasp$wtx_new_inc <- slag(datasp$new_inc, w, 1)
datasp$wtx_pri_en <- slag(datasp$pri_en, w, 1)
datasp$wtx_high_stu <- slag(datasp$high_stu, w, 1)
#模型
## 1.SLX模型(Spatially Lagged X)
## SLX+空间固定效应
FEindivSLX = plm(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu +
                         wtx_kj + wtx_l + wtx_ks + wtx_pe + wtx_inex + wtx_new_inc + wtx_pri_en
                 + wtx_high_stu, data=datasp, model="within", effect="individual")
summary(FEindivSLX)
### 检验空间固定
summary(fixef(FEindivSLX, effect= "individual"))
## SLX+时空固定
FEbothSLX = plm(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu +
                        wtx_kj + wtx_l + wtx_ks + wtx_pe + wtx_inex + wtx_new_inc + wtx_pri_en
                + wtx_high_stu, data=datasp, model="within", effect="twoways")
summary(FEbothSLX)
### 检验时空固定
summary(fixef(FEbothSLX, effect= "individual"))
summary(fixef(FEbothSLX, effect= "time"))

# 2.SLM 模型(可效应分解)
## SLM + 空间固定效应
FEindivslag = spml(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu, data=datasp, listw=w,
                   model="within", effect="individual", lag=TRUE, spatial.error="none")
summary(FEindivslag)
##效应分解
 set.seed(12345)
 library(spatialreg)
 imslm <- impacts(FEindivslag, listw=w, time = 1000)
 summary(imslm, zstats=TRUE, short=F)
## SLM 模型 + 时空双固定效应
 FEslag = spml(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu, data=datasp,
               listw=w, model="within", effect="twoways", lag=TRUE, spatial.error="none")
 summary(FEslag)
# 3.SEM 模型
 ## SEM + 空间固定效应
 FEindivserr = spml(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu, data=datasp,
                    listw=w, model="within", effect="individual", lag=FALSE, spatial.error="kkp")
 summary(FEindivserr)
 ## SEM + 时空双固定效应
 FEserr = spml(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en +
                       high_stu, data = datasp, listw = w, model="within", effect="twoways", lag=FALSE, spatial.error="kkp")
 summary(FEserr)
# 4.SDM 模型
 ## SDM + 空间固定效应
 FEindivDurbin = spml(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en +
                              high_stu, data = datasp, listw = w, model="within", effect="individual", lag=TRUE, spatial.error="none")
 summary(FEindivDurbin)
 ## SDM + 时空双固定效应
 FEDurbin = spml(gdp ~  kj + l + ks + pe + inex + new_inc + pri_en + high_stu +
                         wtx_kj + wtx_l + wtx_ks + wtx_pe + wtx_inex + wtx_new_inc + wtx_pri_en
                 + wtx_high_stu, data=datasp, listw=w, model="within", effect="twoways", lag=TRUE, spatial.error="none")
 summary(FEDurbin)
 
# 六、边际效用
 # SLM估计系数
 library(spatialreg)
 slm <- lagsarlm(gdp2017 ~  kj2017 + l2017 + ks2017 + pe2017 + inex2017 + new_inc2017 + pri_en2017 + high_stu2017,
                 data = cdatashpt, w, method = "eigen", type="lag")
 summary(slm)
 
 # 边际效用计算
 ## 合成自变量X
 X <- cbind(1,cdatashpt$kj2017,cdatashpt$l2017,cdatashpt$ks2017,
            cdatashpt$pe2017,cdatashpt$inex2017,cdatashpt$new_inc2017,cdatashpt$pri_en2017,cdatashpt$high_stu2017)
 head(X)
 ## 变量增加前因变量（gdp2017）的预测值
 # The pre-predicted values
 rho <- slm$rho # Estimated rho from SLM model
 beta_hat <- coef(slm)[-1] # Estimated parameters
 A <- invIrW(w, rho = rho) # (I - rho*W)^{-1}
 X <- cbind(1,cdatashpt$kj2017,cdatashpt$l2017,cdatashpt$ks2017,
            cdatashpt$pe2017,cdatashpt$inex2017,cdatashpt$new_inc2017,
            cdatashpt$pri_en2017,cdatashpt$high_stu2017) # Matrix of observed variables
 y_hat_pre <- A %*% crossprod(t(X), beta_hat) # y hat
 ## 变量增加后因变量（gdp2017）的预测值
 # The post-predicted values
 col_new <- cdatashpt # copy the data frame
 # change the income value
 col_new@data[col_new@data$city == 1, "kj2017"] <- 2.2454180
 # The post-predicted values
 X_d <- cbind(1, col_new$kj2017, cdatashpt$l2017,cdatashpt$ks2017,
              cdatashpt$pe2017,cdatashpt$inex2017,cdatashpt$new_inc2017,cdatashpt$pri_en2017,cdatashpt$high_stu2017)
 y_hat_post <- A %*% crossprod(t(X_d), beta_hat)
 ## 增加前后因变量预测值的差值，即边际效应
 # The difference
 delta_y <- y_hat_post - y_hat_pre
 col_new$delta_y <- delta_y
 # Show the effects
 summary(delta_y)
 sum(delta_y)
 # 可视化边际效应
 library(RColorBrewer)
 library(classInt)
 # 设置颜色主题
 pal6 <- brewer.pal(6, "Oranges")
 cats5 <- classIntervals(col_new$delta_y, n = 5, style = "fisher")
 colors6 <- findColours(cats5, pal6)
 plot(col_new, col = colors6)
 legend("topleft", legend = round(cats5$brks, 5), fill = pal6, bty = "n")
 spplot(col_new, c("delta_y"))
 
 # 七、效应分解
 set.seed(1)
 imslm <- impacts(slm, listw=w, R = 200)
 summary(imslm, zstats = TRUE, short = TRUE)
 
 # 读取MATLAB
 install.packages("R.matlab")
 library(R.matlab) 

 #水刊三板斧
 # 1.导入并转换权重矩阵
 w1 <- read_xlsx("./spatial_analysis/2020fall/data/w1.xlsx",col_names = FALSE)
 w2 <- as.matrix(w1)
 w <- mat2listw(w2,style="W")
 # 2.LM test:空间误差or空间滞后
 library(splm)
 slmtest(FEeffect, listw=w, test="lml")
 slmtest(FEeffect, listw=w, test="rlml")
 slmtest(FEeffect, listw=w, test="lme")
 slmtest(FEeffect, listw=w, test="rlme")
 # 3. SLM 模型 + 空间固定效应(可效应分解)
 # 空间面板滞后项
 datasp$wtx_gdp <- slag(datasp$gdp, w, 1)
 datasp$wtx_kj <- slag(datasp$kj, w, 1)
 datasp$wtx_l <- slag(datasp$l, w, 1)
 datasp$wtx_ks <- slag(datasp$ks, w, 1)
 datasp$wtx_pe <- slag(datasp$pe, w, 1)
 datasp$wtx_inex <- slag(datasp$inex, w, 1)
 datasp$wtx_new_inc <- slag(datasp$new_inc, w, 1)
 datasp$wtx_pri_en <- slag(datasp$pri_en, w, 1)
 datasp$wtx_high_stu <- slag(datasp$high_stu, w, 1)
 
 ## SLX+空间固定效应
 FEindivSLX = plm(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu +
                    wtx_kj + wtx_l + wtx_ks + wtx_pe + wtx_inex + wtx_new_inc + wtx_pri_en
                  + wtx_high_stu, data=datasp, model="within", effect="individual")
 summary(FEindivSLX)
 ### 检验空间固定
 summary(fixef(FEindivSLX, effect= "individual"))
 ## SLX+时空固定
 FEbothSLX = plm(gdp ~ kj + l + ks + pe + inex + new_inc + pri_en + high_stu +
                   wtx_kj + wtx_l + wtx_ks + wtx_pe + wtx_inex + wtx_new_inc + wtx_pri_en
                 + wtx_high_stu, data=datasp, model="within", effect="twoways")
 summary(FEbothSLX)
 
# surprise 爬取新行政代码
 library(rvest)
 library(tidyverse)
 #读取网页 解析并保存为 html 对象
 html <- read_html('http://www.mca.gov.cn/article/sj/xzqh/2020/2020/202003301019.html')
 #提取表格
 html %>% 
   html_table() -> table
 table %>% 
   as.data.frame() %>% 
   as_tibble() %>% 
   select(X2, X3) %>% 
   slice(-1, -2, -3) %>% 
   purrr::set_names("行政区划代码", "单位名称") %>% 
   dplyr::filter(str_length(行政区划代码) == 6) -> mytable
 
 mytable
# 分离省市县
 mytable %>% 
   dplyr::filter(str_sub(行政区划代码, 5, 6)  == "00")
# 省 后四位行政区划代码都是 0
 mytable %>% 
   dplyr::filter(str_sub(行政区划代码, 3, 6)  == "0000")
 # 区县
 mytable %>% 
   dplyr::filter(str_sub(行政区划代码, 5, 6)  != "00")%>% writexl::write_xlsx("temp.xlsx")

 install.packages("writexl")
library(writexl)
 