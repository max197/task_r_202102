#Elaborado por: Max Gómez
#Taller 3: PARTE B

# initial configuration
if (!require("pacman")) install.packages("pacman") # Isntalar pacman (sino está instalada)
require(pacman) # llamar pacman
p_load(tidyverse,viridis,sf,leaflet,raster,maps,Rfast) # llamar y/o instalar librerias
library(skimr)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8") # Encoding UTF-8


#*********
#*PUNTO 1
#*********

###
#1.1 Importar datos espaciales. 
###

#1.1.1

via = st_read("task_3/data/input/VIAS.shp")
puntos = st_read("task_3/data/input/MGN_URB_TOPONIMIA.shp")

#1.1.2

c_medico = filter(puntos, CSIMBOL %in%  c( "021001", "021002", "021003"))

#1.1.3

c_poblado = readRDS("~/TallerR/task1/task_r_202102/task_3/data/input/c poblado (2017).rds")
#Dejar códigos DANE dentro de rango
c_poblado = filter(c_poblado, cod_dane >= 54001 & cod_dane< 55000)

depto  = readRDS("~/TallerR/task1/task_r_202102/task_3/data/input/dp deptos (2017).rds")
#Dejar polígono de Norte de Santander
depto = filter(depto, name_dpto == "NORTE DE SANTANDER")

mapmuse = readRDS("~/TallerR/task1/task_r_202102/task_3/data/input/victimas_map-muse.rds")

##
#1.2 ATRIBUTOS DE LOS OBJETOS
##

#Además de skim, uso head() para ver las primeras entradas de cada objeto
#y así enterarme cómo se ven las entradas de cada columna.
#También uso names() para ver el nombre de las columnas. 

skim(depto)
head(depto) 
#Simple feature collection with 1 feature and 2 fields
#Geometry type: MULTIPOLYGON
#Dimension:     XY
#Bounding box:  xmin: -73.63379 ymin: 6.872201 xmax: -72.04761 ymax: 9.290847

skim(c_poblado)
head(c_poblado)
#Geodetic CRS:  WGS 84

#skim(mapmuse) lo comento porque se demora mucho
names(mapmuse)
table(mapmuse$genero)
head(factor(mapmuse$condicion))
table(mapmuse$condicion)
head(factor(mapmuse$actividad))
table(mapmuse$actividad)

#Si corro skim(via) se estalla el programa. 
#Si corro skim(puntos) se estalla el programa



###
#1.3 GEOMETRÍAS DEL OBJETO
###

#1.3.1
head(mapmuse) #-> Geodetic CRS: WGS 84
st_bbox(mapmuse)

head(c_poblado)  #-> Geodetic CRS: WGS 84
st_bbox(c_poblado) 

head(depto)#-> Geodetic CRS: WGS 84
st_bbox(depto)
 
head(puntos) #-> Geodetic CRS: WGS 84
st_bbox(puntos)

head(c_medico) #-> Geodetic CRS: WGS 84
st_bbox(c_medico)

head(c_poblado) #-> Geodetic CRS: WGS 84
st_bbox(c_poblado)

##1.3.2

#el crs que queremos agregar:
crs ='+proj=utm +zone=19 +datum=WGS84 +units=m +no_def'

#usamos la función st_transform()

c_medico = st_transform(c_medico,crs)
c_poblado = st_transform(c_poblado,crs)
depto = st_transform(depto,crs)
mapmuse = st_transform(mapmuse,crs)
puntos = st_transform(puntos,crs)
via = st_transform(via,crs)

###
#1.4OPERACIONES GEOMÉTRICAS
###

#1.4.1

#clipping: Dejamos los accidentes que ocurrieron en el departamento 'Norte de Santander'
accidentes_santader = st_crop(mapmuse,depto)
plot(accidentes_santader)

#1.4.2

poblado = c_poblado[1,]
vias_poblado = st_crop(via,poblado)

vias_poblado$largo = lapply(vias_poblado$geometry,st_length) #distancia en metros




###
#PUNTO 2
###


#******#******#******#******#******#******
#2.0
#******#******#******#******#******#******

#Crear variables dist_vias, dist_cmedico y dist_cpoblado

#Idea tomada de https://gis.stackexchange.com/questions/351085/r-find-n-closest-points-to-each-point-in-spatialpointsdataframe

dist_vias = st_distance(mapmuse$geometry,via$geometry) #calcula una matriz de distancias
#donde las filas representan los puntos en mapmuse y las columnas en el objeto via.

#rowMins viene del paquete Rfast

mapmuse$dist_vias = rowMins(dist_vias) #Por cada fila, encuentre la distancia del elemento más cercano
#i.e., mire entre las distintas columnas (pero en la misma fila), cual es el el elemento en 'via'
#más cercano al punto en map muse. Esto lo hace rowMins para cada fila. 

#La misma lógica para dist_cmedico y dist_cpoblado:

dist_cmedico = st_distance(mapmuse$geometry,c_medico$geometry)
mapmuse$dist_cmedico = rowMins(dist_cmedico)

dist_cpoblado = st_distance(mapmuse$geometry,c_poblado$geometry)
mapmuse$dist_cpoblado = rowMins(dist_cpoblado)


#Construir la variable fallecido
mapmuse = mutate(mapmuse, fallecido = ifelse(estado=="Muerto",1,0))

p_load(sfheaders)
mapmuse_df = sf_to_df(mapmuse,fill=TRUE)




#******#******#******#******#******#******
#2.1 Estimar un modelo de probabilidad lineal
#******#******#******#******#******#******



p_load(broom, # tidy-coefficients
       mfx, # marginal effects
       margins,  # marginal effects
       estimatr, # robust standard errors
       lmtest, # HAC (Newey-West) standard errors
       fixest, # hdfe regressions (feols)
       modelsummary, # Coefplot with modelplot
       stargazer, # export tables to latex 
       rockchalk #para función outreg
)  

ols = lm(fallecido ~ tipo_accidente + year + condicion + genero + actividad + cod_mpio + dist_vias + dist_cmedico + dist_cpoblado , data = mapmuse) 
ols %>% summary() 

#******#******#******#******#******#******
#2.2 coefplot de Modelo de probabilidad lineal
#******#******#******#******#******#******

coefplot(ols)


#******#******#******#******#******#******
#2.3 Estimar un modelo de probit y uno logit
#******#******#******#******#******#******

logit = glm(fallecido ~ tipo_accidente + year + condicion + genero + actividad + cod_mpio + dist_vias + dist_cmedico + dist_cpoblado , data = mapmuse , family = binomial(link="logit")) 
logit %>% summary()

probit = glm(fallecido ~ tipo_accidente + year + condicion + genero + actividad + cod_mpio + dist_vias + dist_cmedico + dist_cpoblado , data = mapmuse , family = binomial(link="probit")) 
probit %>% summary()


#******#******#******#******#******#******
#2.4
#******#******#******#******#******#******
outreg(list("Linear PM" = ols, "Probit" = probit, "Logit" =logit), type="html",browse=TRUE) 

#******#******#******#******#******#******
#2.5
#******#******#******#******#******#******

p_load(sjPlot) #para usar función plot_model que pinta los efectos marginales

plot_model(probit, type="pred", term= c('dist_cmedico')) #pinta los marginal effects
# de la distancia al centro médicosobre sobre la prob. de fallecer
plot_model(logit, type="pred", term= c('dist_cmedico'))
           



###
#PUNTO 3
###


#******#
#3.1
#******#

library(xml2) #para leer html
library(rvest)

url = "https://es.wikipedia.org/wiki/Departamentos_de_Colombia"
wikipedia = read_html(url)

#******#
#3.2
#******#

titulo_pagina = wikipedia %>% html_nodes(xpath = '//*[@id="firstHeading"]') %>% html_text()


#******#
#3.3
#******#

tabla = wikipedia %>% html_nodes(xpath = '/html/body/div[3]/div[3]/div[5]/div[1]/table[3]') %>% html_table()

