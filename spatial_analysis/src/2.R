#  Learn R with Dr.Hu and His Friends
#  Sun Zhaoyang R lecture No.2 2022.03.11
#  2022 Spring

# 请点击以下两个链接、注册账号
# 自然资源部 http://bzdt.ch.mnr.gov.cn/
# 高德地图API https://console.amap.com/dev/id/phone


#  set current folder as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#  setwd("D://xxx//xxx")

# 国家自然资源部 标准地图服务系统
# Chrome浏览器打开
# http://bzdt.ch.mnr.gov.cn/

library(sf)
library(ggplot2)


################
# 0. 画中国地图的常见问题
################
# 1 中国台湾部分的缺失
# 2 南海与九段线的缺失
# 3 西藏交界的中印边界划分有误
# 4 新疆与西藏交界的中印边界划分有误


# 中华人民共和国民政部 全国行政区划信息查询平台
# 提供了省级与县级两种类型的地图，其审图号为：GS(2018)2512 
# 可以查看民政部官网的源代码，点击请求，返回json格式地图数据
# API前缀是 http://xzqh.mca.gov.cn/data/,
# 审图号：GS(2018)2512号

# 获取全国省级地图，则加后缀quanguo.json
# 获取全国县级地图，则加后缀xian_quanguo.json
# 获取部分地区，如某个市的县级地图，则加该行政区域代码，再加.json
# 如果要获取市级地图，需要按遍历行政区域代码获取所有市的地图，然后合并县级区域
# 全国主要山脉，南海九段线数据，则加后缀quanguo_Line.geojson

API_pre = "http://xzqh.mca.gov.cn/data/"


################
# 1. 全国地图
################
China <- st_read(dsn = paste0(API_pre, "quanguo.json"),stringsAsFactors=FALSE)
China_county <- st_read(dsn = paste0(API_pre, "xian_quanguo.json"),stringsAsFactors=FALSE)

ggplot()+ # 底图
  geom_sf(data = China, aes(fill = QUHUADAIMA)) # 省域

ggplot()+ # 底图
  geom_sf(data = China, aes(fill = FillColor)) # 省域

ggplot()+ # 底图
  geom_sf(data = China_county) # 县域
ggplot()+ # 底图
  geom_sf(data = China_county, aes(fill = FillColor)) # 县域

ggplot()+ # 底图
  geom_sf(data = China, aes(fill = QUHUADAIMA))+ # 省域
  coord_sf()  # 投影方式

# 全国法院
library(readxl)
dat_court <- read_excel("D:\\AAA所有政治学研究\\2022律所与法院\\全国法院信息 审理法院字典 20191210.xlsx", sheet = "新数据源")
dat_lawfirm <- read_excel("D:\\AAA所有政治学研究\\2022律所与法院\\法律相关poi.xlsx", sheet="法律相关poi")

ggplot()+ # 底图
  geom_sf(data = China, alpha = 0.5)+ # 省域
  geom_point(data = dat_court, aes(x = 经度, y = 纬度), size=0.001, alpha = 0.8, color = 'red')+
  geom_point(data = dat_lawfirm, aes(x = display_x, y = display_y), size=0.001, alpha = 0.8, color = 'grey')+
  coord_sf()+
  theme_bw()


################
# 2. 国境线
################
China_line = st_read(dsn = paste0(API_pre, "quanguo_Line.geojson"), 
                     stringsAsFactors=FALSE) 
Border_line <- China_line[China_line$QUHUADAIMA == "guojiexian",]

JiuDuanXian <- subset.data.frame(Border_line,)

ggplot()+ # 底图
  geom_sf(data = Border_line)

################
# 3. 自主调色
################
ColorManual <-as.numeric(China$NAME)

zhuose_data <- read.csv("https://raw.githubusercontent.com/slyang-cn/data/slyangcn/your_data.csv")
zhuose_data$QUHUADAIMA <- as.character(zhuose_data$QUHUADAIMA) # 因China数据中QUHUADAIMA是chr类型
CHINA <- dplyr::left_join(China, zhuose_data, by= "QUHUADAIMA") 

ggplot()+ 
  geom_sf(data = CHINA, aes(fill = yanse))

ggplot()+
  geom_sf(data = CHINA, aes(fill = factor(yanse))) +
  scale_fill_manual("class", values=c("#FFCCCC", "#FF9333", "#FF6660","#FF5111","#CC0070"),
                    breaks = c("0~200","200~400","400~600","600~1000","1000+"),
                    labels = c("0~200","200~400","400~600","600~1000","1000+"))


################
# 4. 世界地图
################
library(mapdata)

# 绘制基本地图
mapworld<-borders("world", colour = "gray50", fill = "white") 
mapworld

ggplot()+
  mapworld+xlim(-180,0)+ylim(-90,90) # 利用ggplot呈现，同时地图纵坐标范围从-60到90


# 实用技术 | 如何用R绘制并填充相对正确的世界地图
# https://mp.weixin.qq.com/s/m-jF0P2jGm7o6GQ3-YLBlQ

# 最省事的方法

# 国家自然资源部 标准地图服务系统 Chrome浏览器打开
# http://bzdt.ch.mnr.gov.cn/


################
# 5. 网络地图信息获取
################
library(tidyverse)
library(magrittr)
library(jsonlite)

# 对fromJSON()函数进行改造，使之仅返回经纬度信息
lonlat <- function(x) return(fromJSON(x)$geocodes$location)

# 示例地点
data <- data.frame(dizhi = c("清华大学", "展览路加油站", "什刹海"))

## 获取地理信息
key <- "c5335ae546059d163521c1f1ac71e60c"  ## 替换成自己的Key
data %<>%
  mutate(url = paste0("https://restapi.amap.com/v3/geocode/geo?",
                      "key=", key,
                      "&address=", dizhi),json = map(url, lonlat))
data$json

# 详细教程
# https://blog.csdn.net/weixin_54000907/article/details/123321115


################
# 6. 直接调用网络地图
################
library(leaflet)
leaflet()
names(providers)# 地图服务提供商

# 添加瓦片地图
leaflet() %>%
  addTiles()# 默认的OpenStreetMap地图

# 加载局部地图
setView(map, lng, lat,
        zoom, options = list())

# 设置观察范围
leaflet() %>%
  addTiles() %>%
  setView(lng = 116.35, lat = 39.97, zoom = 12)

leaflet() %>%
  addTiles() %>%
  fitBounds(-74,39,-75,40)

# 添加动画效果
leaflet() %>%
  addTiles() %>%
  flyTo(lng = -75, lat = 39.9, zoom = 10)

# 使用中国出品的地图
library(leafletCN)
leaflet() %>%
  setView(lng = 116.35, lat = 39.97, zoom = 12) %>%
  amap() # 比如高德


# 往地图对象中添加一个标志物
library(leafletCN)
m <- leaflet() %>%
  amap() %>%
  setView(lng = 116.35, lat = 39.97, zoom = 12) 
m %>%
  addMarkers(lng = 116.344639, lat = 39.930913)

m %>%
  addPopups(lng = 116.344639, lat = 39.930913,
            popup = "展览路加油站")

# 自主设置标记点
doge <- makeIcon(iconUrl = "doge.png",
                 iconWidth = 40, iconHeight = 40)  

# 批量添加标记点：第一种方法
m %>%
  addMarkers(lng = c(116.339996, 116.344639, 116.392555),
             lat = c(39.985401, 39.930913, 39.933781),
             icon = doge)   

# 批量添加标记点：第二种方法
dta <- data.frame(x = c(116.339996, 116.344639, 116.392555),
                  y = c(39.985401, 39.930913, 39.933781))

m %>%
  addMarkers(lng = ~x, lat = ~y,
             data = dta,
             icon = doge)  

# 添加图形
m %>%
  addCircles(
    map,
    lng = 116.344639,
    lat = 39.930913,
    radius = 3000,
    group = NULL,
    stroke = TRUE,
    color = "grey",
    weight = 5,
    opacity = 0.5,
    fill = TRUE,
    fillColor = "grey",
    fillOpacity = 0.05) %>%
  addPopups(
    lng = 116.344639, lat = 39.930913, popup = "展览路")

################
# Appendix. 一些可能有用的链接
################

# 添加特殊符号文本：https://tidyfriday.cn/posts/17200/
# 改变坐标轴的设置：https://statisticsglobe.com/change-font-size-of-ggplot2-plot-in-r-axis-text-main-title-legend
# 颜色设置方式其一：https://www.rapidtables.com/web/color/RGB_Color.html
# 颜色设置方式其二：http://sape.inf.usi.ch/quick-reference/ggplot2/colour

