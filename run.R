source('R/fgwc.R')
source('R/ei.R')
source('R/gsafgwc.R')
source('R/ifafgwc.R')
source('R/hhofgwc.R')
source('R/psofgwc.R')
source('R/tlbofgwc.R')
source('R/abcfgwc.R')
source('R/fpafgwc.R')

a = 1
if (a !=1) print ('babi') else print('ayam')
library(roxygen2)
roxygenise()
?fgwcuv
library(naspaclust)
data('census2010')
data('census2010dist')
data('census2010pop')
b = fgwcuv(census2010,census2010pop,census2010dist,'u','3',2,'euclidean',4)

options(warn=-1)
b = fgwc(census2010,census2010pop,census2010dist,'classic',c(1),c(1))
b = fgwc(census2010,census2010pop,census2010dist,'hho',c(1),c(1))

cc = hhofgwc(census2010,census2010pop,census2010dist)


param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
              alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
## tune the FPA parameter
fpa_param <- c(vi.dist='normal',npar=5,same=15,p=0.7,gamma=1.2,lambda=1.5,ei.distr='logchaotic',chaos=3)
##FGWC with FPA
res2 = fgwc(census2010,census2010pop,census2010dist,'fpa',param_fgwc,fpa_param)
