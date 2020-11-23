fgwc <- function(data,pop,distmat,algorithm='classic',fgwc_param,opt_param){
    if (!fgwc_param['kind']%in%c('u','v','NA')) kind <- NA else kind <- fgwc_param['kind']
    if (is.na(fgwc_param['ncluster'])|fgwc_param['ncluster']<2) ncluster <- 2 else ncluster <- as.numeric(fgwc_param['ncluster'])
    if ((is.na(fgwc_param['m'])|fgwc_param['m']<1)) m <- 2 else m <- as.numeric(fgwc_param['m'])
    if ((is.na(fgwc_param['distance']))) distance <- 'euclidean' else distance <- as.numeric(fgwc_param['distance'])
    if ((is.na(fgwc_param['order'])|fgwc_param['order']<=0)) order <- 2 else order <- as.numeric(fgwc_param['order'])
    if ((is.na(fgwc_param['alpha'])|fgwc_param['alpha']>1|fgwc_param['alpha']<0)) alpha <- 0.7 else alpha <- as.numeric(fgwc_param['alpha'])
    if ((is.na(fgwc_param['a'])|fgwc_param['a']<0)) a <- 1 else a <- as.numeric(fgwc_param['a'])
    if ((is.na(fgwc_param['b'])|fgwc_param['b']<0)) b <- 1 else b <- as.numeric(fgwc_param['b'])
    if ((is.na(fgwc_param['max.iter'])|fgwc_param['max.iter']<0)) max.iter <- 500 else max.iter <- as.numeric(fgwc_param['max.iter'])
    if ((is.na(fgwc_param['error'])|fgwc_param['error']<0)) error <- 1e-5 else error <- as.numeric(fgwc_param['error'])
    if ((is.na(fgwc_param['randomN'])|fgwc_param['randomN']<0)) randomN <- 1 else randomN <- as.numeric(fgwc_param['randomN'])
    
    if(algorithm!='classic'){
        if(is.na(opt_param['vi.dist'])) vi.dist <- 'uniform' else vi.dist <- opt_param['vi.dist']
        if(is.na(opt_param['npar'])|opt_param['npar']<0) npar <- 10 else npar <- as.numeric(opt_param['npar'])
        if(is.na(opt_param['par.no'])|opt_param['par.no']<0) par.no <- 2 else par.no <- as.numeric(opt_param['par.no'])
        if(is.na(opt_param['par.dist'])|opt_param['par.dist']<0) par.dist <- 'euclidean' else par.dist <- opt_param['par.dist']
        if(is.na(opt_param['par.order'])|opt_param['par.order']<0) par.order <- 2 else par.order <- as.numeric(opt_param['par.order'])
        if(is.na(opt_param['pso'])|opt_param['pso']<0) pso <- T else pso <- as.logical(opt_param['pso'])
        if(is.na(opt_param['same'])|opt_param['same']<0) same <- 10 else same <- as.numeric(opt_param['same'])
        if(is.na(opt_param['type'])) type <- 'sim.annealing' else type <- opt_param['type']
        if(is.na(opt_param['ei.distr'])) ei.distr <- 'normal' else ei.distr <- opt_param['ei.distr']
        if(is.na(opt_param['vmax'])|opt_param['vmax']<0) vmax <- 0.7 else vmax <- as.numeric(opt_param['vmax'])
        if(is.na(opt_param['wmax'])|opt_param['wmax']<0) wmax <- 0.9 else wmax <- as.numeric(opt_param['wmax'])
        if(is.na(opt_param['wmin'])|opt_param['wmin']<0) wmin <- 0.4 else wmin <- as.numeric(opt_param['wmin'])
        if(is.na(opt_param['chaos'])|opt_param['chaos']<0) chaos <- 4 else chaos <- as.numeric(opt_param['chaos'])
        if(is.na(opt_param['x0'])|opt_param['x0']<0) x0 <- 'F' else x0 <- opt_param['x0']
        if(is.na(opt_param['map'])|opt_param['map']<0) map <- 0.7 else map <- as.numeric(opt_param['map'])
        if(is.na(opt_param['ind'])|opt_param['ind']<0) ind <- 1 else ind <- as.numeric(opt_param['ind'])
        if(is.na(opt_param['skew'])|opt_param['skew']<0) skew <- 0 else skew <- as.numeric(opt_param['skew'])
        if(is.na(opt_param['sca'])|opt_param['sca']<0) sca <- 1 else sca <- as.numeric(opt_param['sca'])
    }

    if (algorithm=='classic'){
        return(fgwcuv(data, pop, distmat, kind=kind,ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, max.iter=max.iter, error=error,randomN=randomN, uij=NA, vi=NA))
    }
    else if(algorithm=='abc'){
        opt_param <- get_param_abc(opt_param)
        return(abcfgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN, vi.dist=vi.dist, 
                      nfood=npar, n.onlooker=as.numeric(opt_param['n.onlooker']), 
                      limit=opt_param['limit'], pso=pso, abc.same=same))
    }
    else if(algorithm=='gsa'){
        opt_param <- get_param_gsa(opt_param)
        return(gsafgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter,randomN=randomN,
                      vi.dist=vi.dist,npar=npar,par.no=par.no,par.dist=par.dist, par.order=par.order,
                      gsa.same=same, G=opt_param['G'], g.type=type,vmax=vmax, pso=pso,wmax=wmax,wmin=wmin,
                      chaos=chaos,x0=x0,map=map,ind=ind,skew=skew,sca=sca))
    }
    else if(algorithm=='fpa'){
        opt_param <- get_param_fpa(opt_param)
        return(fpafgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist, nflow=npar, p=opt_param['p'], gamma=opt_param['gamma'], 
                      lambda=opt_param['lambda'], delta=opt_param['delta'],
                      ei.distr=et.distr,flow.same=same,r=chaos,m.chaotic=map,skew=skew,sca=sca))
    }
    else if(algorithm=='hho'){
        opt_param <- get_param_hho(opt_param)
        return(hhofgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist,nhh=npar,hh.alg=opt_param['algo'],
                      A=as.numeric(opt_param[c('a1','a2','a3')]),p=as.numeric(opt_param['p']),
                      hh.same=same,levy.beta=as.numeric(opt_param['beta']),as.numeric(opt_param['update.type'])))
    }
    else if(algorithm=='ifa'){
        opt_param <- get_param_ifa(opt_param)
        return(ifafgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist, ei.distr=ei.distr,fa.same=same, nfly=npar, ffly.no=par.no, 
                      ffly.dist=par.dist, ffly.order=par.order, gamma=as.numeric(opt_param['gamma']), 
                      ffly.beta=as.numeric(opt_param['beta']),ffly.alpha=as.numeric(opt_param['alpha']),
                      r.chaotic=chaos,m.chaotic=map,ind.levy=ind,skew.levy=skew,
                      scale.levy=sca,ffly.alpha.type=as.numeric(opt_param['update_type'])))
    }
    else if(algorithm=='pso'){
        opt_param <- get_param_pso(opt_param)
        return(psofgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist,npar=npar,vmax=vmax, pso.same=same, c1=as.numeric(opt_param['c1']),
                      c2=as.numeric(opt_param['c2']), w.inert=type,
                      wmax=wmax,wmin=wmin,chaos=chaos,x0=x0,map=map,ind=ind,skew=skew,sca=sca))
    }
    else if(algorithm=='tlbo'){
        opt_param <- get_param_tlbo(opt_param)
        return(tlbofgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist,nstud=npar,vmax=vmax, tlbo.same=same,
                      nselection=as.numeric(opt_param['nselection']),
                      elitism=as.logical(opt_param['elitism']),n.elite=as.numeric(opt_param['n.elite'])))
    }
    
}

get_param_abc <- function(param){
    if(is.na(param['n.onlooker'])|param['n.onlooker']<0) param['n.onlooker'] <- 5
    if(is.na(param['limit'])|param['limit']<0) param['limit'] <- 4
    return (param)
}
get_param_gsa <- function(param){
    if(is.na(param['G'])|param['G']<0) param['G'] <- 1
    return(param)
}

get_param_fpa <- function(param){
    if(is.na(param['p'])|param['p']<0|param['p']>1) param['p'] <- 0.8
    if(is.na(param['gamma'])|param['gamma']<0) param['gamma'] <- 1
    if(is.na(param['lambda'])|param['lambda']<0) param['lambda'] <- 1.5
    if(is.na(param['delta'])|param['delta']<0) param['delta'] <- 0
    return(param)
}

get_param_hho <- function(param){
    if(is.na(param['algo'])|param['algo']<0) param['algo'] <- 'heidari'
    if(is.na(param['a1'])|param['a1']<0) param['a1'] <- 2
    if(is.na(param['a2'])|param['a2']<0) param['a2'] <- 1
    if(is.na(param['a3'])|param['a3']<0) param['a3'] <- 0.5
    if(is.na(param['p'])|param['p']<0|param['p']>1) param['p'] <- 0.5
    if(is.na(param['beta'])|param['beta']<0) param['gamma'] <- 1
    if(is.na(param['update_type'])|param['update_type']<5) param['update_type'] <- 5
    return(param)
}

get_param_ifa <- function(param){
    if(is.na(param['gamma'])|param['gamma']<0) param['gamma'] <- 1
    if(is.na(param['beta'])|param['beta']<0) param['gamma'] <- 1
    if(is.na(param['alpha'])|param['alpha']<0) param['alpha'] <- 1
    if(is.na(param['update_type'])|param['update_type']<5) param['update_type'] <- 4
    return(param)
}

get_param_pso <- function(param){
    if(is.na(param['c1'])|param['c1']<0) param['c1'] <- 0.49
    if(is.na(param['c2'])|param['c2']<0) param['c2'] <- 0.49
    return(param)
}

get_param_tlbo <- function(param){
    if(is.na(param['nselection'])|param['nselection']<0) param['nselection'] <- 10
    if(is.na(param['elitism'])|param['elitism']<0) param['elitism'] <- F
    if(is.na(param['n.elite'])|param['n.elite']<0) param['n.elite'] <- 2
    return(param)
}