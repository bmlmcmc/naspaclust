source('R/fgwc.R')
source('R/mainfunction.R')
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


## tune the ABC parameter
abc_param <- c(vi.dist='normal',npar=10,pso=FALSE,same=15,n.onlooker=5,limit=5)
## FGWC with ABC optimization algorithm
res2 = fgwc(census2010,census2010pop,census2010dist,'abc',param_fgwc,abc_param)

param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
              alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
## tune the FPA parameter
fpa_param <- c(vi.dist='normal',npar=5,same=15,p=0.7,gamma=1.2,lambda=1.5,ei.distr='logchaotic',chaos=3)
##FGWC with FPA

res2 = fgwc(census2010,census2010pop,census2010dist,'fpa',param_fgwc,fpa_param)

gsa_param <- c(vi.dist='normal',npar=5,same=15,G=1,vmax=0.7,new=F)
##FGWC with GSA
res2 = fgwc(census2010,census2010pop,census2010dist,'gsa',param_fgwc,gsa_param)


## tune the HHO parameter
hho_param <- c(vi.dist='normal',npar=5,same=15,algo='bairathi',a1=3,a2=1,a3=0.4)
##FGWC with HHO
res2 = fgwc(census2010,census2010pop,census2010dist,'hho',param_fgwc,hho_param)

ifa_param <- c(vi.dist='uniform', ei.distr='logchaotic',
               fa.same=10, nfly=15, ffly.no=3, ffly.dist='minkowski', ffly.order=4, gamma=1, ffly.beta=1.5,
               ffly.alpha=1, r.chaotic=4,update_type=4)
##FGWC with IFA
source('R/mainfunction.R')
res2 = fgwc(census2010,census2010pop,census2010dist,'ifa',param_fgwc,gsa_param)

## tune the PSO parameter
pso_param <- c(vi.dist='uniform',npar=15,
               vmax=0.8, pso.same=10, c1=0.7, c2=0.6, w.inert='chaotic',
               wmax=0.8,wmin=0.3)
##FGWC with PSO
res2 = fgwc(census2010,census2010pop,census2010dist,'gsa',param_fgwc,pso_param)

## tune the TLBO parameter
tlbo_param <- c(vi.dist="normal",nstud=10,vmax=0.4, tlbo.same=10,
                nselection=10,elitism=F,n.elite=2)
##FGWC with TLBO
res2 = fgwc(census2010,census2010pop,census2010dist,'gsa',param_fgwc,tlbo_param)



a0 = fgwcuv(census2010,census2010pop,census2010dist,'v',3,2,'euclidean',alpha=1)
a1 = fgwcuv(census2010,census2010pop,census2010dist,'v',3,2,'euclidean',)
a2 = gsafgwc(census2010,census2010pop,census2010dist,3,2,'euclidean')
a3 = ifafgwc(census2010,census2010pop,census2010dist,3,2,'euclidean')
a4 = hhofgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',hh.alg = 'bairathi')
a5 = hhofgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',hh.alg = 'heidari')
a6 = psofgwc(census2010,census2010pop,census2010dist,3,2,'euclidean')
a7 = tlbofgwc(census2010,census2010pop,census2010dist,3,2,'euclidean')
a8 = tlbofgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',elitism = T,n.elite = 1)
a9 = abcfgwc(census2010,census2010pop,census2010dist,3,2,'euclidean')
a10 = abcfgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',pso=T)
a11 = fpafgwc(census2010,census2010pop,census2010dist,3,2,'euclidean')

file.edit(".Rprofile")
library(devtools)
use_gpl3_license(version = 3, include_future = TRUE)
build_manual()

Sys.setenv(PATH=paste(Sys.getenv("PATH"),"C:/Users/dokuganryu/AppData/Local/Programs/MiKTeX/miktex/bin/",sep=";"))
