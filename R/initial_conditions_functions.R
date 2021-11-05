add = function(x,t,n,d){
    idx = 1:min(length(d),length(x)-t)
    x[t+idx-1] = x[t+idx-1] + n*d[idx]
    return(x)
}

solve_diff_eqns = function(exposures,N,nt,q,r,s,nS_V,y,dE,dIp,dIs,dIa){
     
     na = length(N)
     
     S = matrix(nrow = na,ncol = nt+1)
     V = matrix(nrow = na,ncol = nt+1)
     E = matrix(nrow = na,ncol = nt+1)
     L = matrix(nrow = na,ncol = nt+1)
     Ia = matrix(nrow = na,ncol = nt+1)
     Ip = matrix(nrow = na,ncol = nt+1)
     Is = matrix(nrow = na,ncol = nt+1)
     R = matrix(nrow = na,ncol = nt+1)
     
     # Set initial conditions
     S[,1] = N
     V[,1] = rep(0,na)
     E[,1] = rep(0,na)
     L[,1] = rep(0,na)
     Ia[,1] = rep(0,na)
     Ip[,1] = rep(0,na)
     Is[,1] = rep(0,na)
     R[,1] = rep(0,na)
     
     # Convolve exposures to get numbers entering each state
     nS_E = matrix(0,nrow = na,ncol = nt)
     nV_E = matrix(0,nrow = na,ncol = nt)
     nV_L = matrix(0,nrow = na,ncol = nt)
     nV_R = matrix(0,nrow = na,ncol = nt)
     nE_I = matrix(0,nrow = na,ncol = nt)
     nE_Ip = matrix(0,nrow = na,ncol = nt)
     nE_Ia = matrix(0,nrow = na,ncol = nt)
     nL_Ia = matrix(0,nrow = na,ncol = nt)
     nIp_Is = matrix(0,nrow = na,ncol = nt)
     nIs_R = matrix(0,nrow = na,ncol = nt)
     nIa_R = matrix(0,nrow = na,ncol = nt)
     # for (a in 1:na){
     #     nE_I[a,] = disc_conv(exposures[a,],dE)
     #     nE_Ip[a,] = y[a]*nE_I[a,]
     #     nE_Ia[a,] = (1-y[a])*nE_I[a,]
     #     nL_Ia[a,] = disc_conv()
     #     nIp_Is[a,] = disc_conv(nE_Ip[a,],dIp)
     #     nIs_R[a,] = disc_conv(nIp_Is[a,],dIs)
     #     nIa_R[a,] = disc_conv(nE_Ia[a,],dIa)
     # }
     
     # Iterate system to update numbers in each state
     for (t in 1:nt){
         for (a in 1:na){
             if (S[a,t]==0 & V[a,t]==0){
                 nS_E[a,t] = 0
                 nV_E[a,t] = 0
                 nV_L[a,t] = 0
             } else {
                 nS_E[a,t] = S[a,t]/(S[a,t]+(q[a,t]+r[a,t]+s[a,t])*V[a,t])*exposures[a,t]
                 nV_E[a,t] = q[a,t]*V[a,t]/(S[a,t]+(q[a,t]+r[a,t]+s[a,t])*V[a,t])*exposures[a,t]
                 nV_L[a,t] = r[a,t]*V[a,t]/(S[a,t]+(q[a,t]+r[a,t]+s[a,t])*V[a,t])*exposures[a,t]
                 nV_R[a,t] = s[a,t]*V[a,t]/(S[a,t]+(q[a,t]+r[a,t]+s[a,t]*V[a,t]))*exposures[a,t]
             }
             nE_I[a,] = add(nE_I[a,],t,nS_E[a,t]+nV_E[a,t],dE)
             nE_Ip[a,] = y[a]*nE_I[a,]
             nE_Ia[a,] = (1-y[a])*nE_I[a,]
             nL_Ia[a,] = add(nL_Ia[a,],t,nV_L[a,t],dE)
             nIp_Is[a,] = add(nIp_Is[a,],t,nE_Ip[a,t],dIp)
             nIs_R[a,] = add(nIs_R[a,],t,nIp_Is[a,t],dIs)
             nIa_R[a,] = add(nIa_R[a,],t,nE_Ia[a,t],dIa)
             
             # S[a,t+1] = max(S[a,t] - nS_E[a,t] - nS_V[a,t],0) # + nV_S[a,t] + nR_S[a,t]
             # V[a,t+1] = V[a,t] + nS_V[a,t] - nV_E[a,t] - nV_L[a,t] #- nV_R[a,t] - nV_S[a,t]
             # E[a,t+1] = E[a,t] + nS_E[a,t] + nV_E[a,t] - nE_Ip[a,t] - nE_Ia[a,t]
             # L[a,t+1] = L[a,t] + nV_L[a,t] - nL_Ia[a,t]
             # Ia[a,t+1] = Ia[a,t] + nE_Ia[a,t] + nL_Ia[a,t] - nIa_R[a,t]
             # Ip[a,t+1] = Ip[a,t] + nE_Ip[a,t] - nIp_Is[a,t]
             # Is[a,t+1] = Is[a,t] + nIp_Is[a,t] - nIs_R[a,t]
             # R[a,t+1] = R[a,t] + nIa_R[a,t] + nIs_R[a,t] #+ nV_R[a,t] - nR_S[a,t]             
         }
         
         S[,t+1] = pmax(S[,t] - nS_E[,t] - nS_V[,t],0) # + nV_S[,t] + nR_S[,t]
         V[,t+1] = V[,t] + nS_V[,t] - nV_E[,t] - nV_L[,t] - nV_R[,t] #- nV_S[,t]
         E[,t+1] = E[,t] + nS_E[,t] + nV_E[,t] - nE_Ip[,t] - nE_Ia[,t]
         L[,t+1] = L[,t] + nV_L[,t] - nL_Ia[,t]
         Ia[,t+1] = Ia[,t] + nE_Ia[,t] + nL_Ia[,t] - nIa_R[,t]
         Ip[,t+1] = Ip[,t] + nE_Ip[,t] - nIp_Is[,t]
         Is[,t+1] = Is[,t] + nIp_Is[,t] - nIs_R[,t]
         R[,t+1] = R[,t] + nIa_R[,t] + nIs_R[,t] + nV_R[,t] #- nR_S[,t]   
     }
     
     return(list(S=S[,1:nt],V=V[,1:nt],E=E[,1:nt],L=L[,1:nt],Ia=Ia[,1:nt],
                 Ip=Ip[,1:nt],Is=Is[,1:nt],R=R[,1:nt],nS_E=nS_E,nS_V=nS_V,
                 nV_E=nV_E,nV_L=nV_L,nV_R=nV_R,nE_I=nE_I,nE_Ip=nE_Ip,nE_Ia=nE_Ia,
                 nL_Ia=nL_Ia,nIp_Is=nIp_Is,nIs_R=nIs_R,nIa_R=nIa_R))
}

calc_init_condns = function(inc_dt,vax_dt_wide,agegroups_model,covy,vrnt_prop,ve_params,dE,dIp,dIs,dIa){
    inc_dt1 = copy(inc_dt)
    
    countries = inc_dt1[,unique(country)]
    ncountries = length(countries)
    
    # Cast variant proportion data table to wide format
    vrnt_prop_wide = dcast(vrnt_prop,country+date~vrnt,value.var = "prop_vrnt")
    cols = paste0("prop_vrnt",c("","2","3"))
    setnames(vrnt_prop_wide,c("Other","Alpha","Delta"),cols)
    
    # Merge with overall data table
    inc_dt1 = merge(inc_dt1,vrnt_prop_wide,by=c("country","date"),all.x=T)
    # Backfill variant proportions with first non-NA observation
    inc_dt1[,(cols):=nafill(.SD,type="nocb"),.SDcols=cols,by=.(country)]
    
    # Merge with vaccination data
    inc_dt1 = merge(inc_dt1,vax_dt_wide[,!"population"],by=c("country","age_group_model","date"),all.x=T)
    
    # Calculate relative proportions of each vaccine type
    inc_dt1[,`:=`(p_va=prop_va/(prop_va+prop_vb),p_vb=prop_vb/(prop_va+prop_vb))]
    
    # Calculate average vaccine efficacy over time against 3 variants, 
    # accounting for changing vaccine type proportions and variant proportions
    inc_dt1[,`:=`(ei=fifelse(!(prop_va==0 & prop_vb==0),
                            p_va*(prop_vrnt*ve_params$ei_va2+prop_vrnt2*ve_params$ei2_va2+prop_vrnt3*ve_params$ei3_va2)+
                                p_vb*(prop_vrnt*ve_params$ei_vb2+prop_vrnt2*ve_params$ei2_vb2+prop_vrnt3*ve_params$ei3_vb2),0),
                 ed=fifelse(!(prop_va==0 & prop_vb==0),
                            p_va*(prop_vrnt*ve_params$ed_va2i+prop_vrnt2*ve_params$ed_va2i2+prop_vrnt3*ve_params$ed_va2i3)+
                                p_vb*(prop_vrnt*ve_params$ed_vb2i+prop_vrnt2*ve_params$ed_vb2i2+prop_vrnt3*ve_params$ed_vb2i3),0),
                 eh=fifelse(!(prop_va==0 & prop_vb==0),
                            p_va*(prop_vrnt*ve_params$eh_va2d+prop_vrnt2*ve_params$eh_va2d2+prop_vrnt3*ve_params$eh_va2d3)+
                                p_vb*(prop_vrnt*ve_params$eh_vb2d+prop_vrnt2*ve_params$eh_vb2d2+prop_vrnt3*ve_params$eh_vb2d3),0),
                 em=fifelse(!(prop_va==0 & prop_vb==0),
                            p_va*(prop_vrnt*ve_params$em_va2d+prop_vrnt2*ve_params$em_va2d2+prop_vrnt3*ve_params$em_va2d3)+
                                p_vb*(prop_vrnt*ve_params$em_vb2d+prop_vrnt2*ve_params$em_vb2d2+prop_vrnt3*ve_params$em_vb2d3),0),
                 et=fifelse(!(prop_va==0 & prop_vb==0),
                            p_va*(prop_vrnt*ve_params$et_va2i+prop_vrnt2*ve_params$et_va2i2+prop_vrnt3*ve_params$et_va2i3)+
                                p_vb*(prop_vrnt*ve_params$et_vb2i+prop_vrnt2*ve_params$et_vb2i2+prop_vrnt3*ve_params$et_vb2i3),0))]
    inc_dt1[,`:=`(q=(1-ei)*(1-ed)*(1-et),r=(1-ei)*ed*(1-et),s=(1-ei)*et)]
    
    prev_list = vector("list",ncountries)
    for (i in 1:ncountries){
        cntry = countries[i]
        exposures = as.matrix(dcast(inc_dt1[country==cntry],age_group_model ~ date,value.var = "exposures")[,-1])
        population = unique(inc_dt1[country==cntry,.(age_group_model,population)])[,population]
        nt = ncol(exposures)
        nS_V = as.matrix(dcast(vax_dt_wide[country==cntry,.(date,age_group_model,nS_V=nS_Va1+nS_Vb1)],age_group_model ~ date,value.var = "nS_V")[,-1])[,1:nt]
        q = as.matrix(dcast(inc_dt1[country==cntry],age_group_model ~ date,value.var = "q")[,-1])
        r = as.matrix(dcast(inc_dt1[country==cntry],age_group_model ~ date,value.var = "r")[,-1])
        s = as.matrix(dcast(inc_dt1[country==cntry],age_group_model ~ date,value.var = "s")[,-1])
        
        Y = solve_diff_eqns(exposures,population,nt,q,r,s,nS_V,covy,dE,dIp,dIs,dIa)
        
        age_grp = agegroups_model[reshape2::melt(Y[[1]])$Var1]
        dts = inc_dt1[country==cntry,unique(date)][reshape2::melt(Y[[1]])$Var2] 
        for (j in seq_along(Y)){
            Y[[j]] = reshape2::melt(Y[[j]],value.name = names(Y)[j])[names(Y)[j]]
        }
        Y = do.call(cbind,Y)
        setDT(Y)
        Y[,`:=`(country=cntry,age_group_model=age_grp,date=dts)]
        
        prev_list[[i]] = Y
    }
    prev_dt = rbindlist(prev_list)
    prev_dt[,age_group_model:=factor(age_group_model,levels=agegroups_model)]
    
    prev_dt = merge(inc_dt1,prev_dt,by=c("country","age_group_model","date"))
    
    return(prev_dt)
}

rep_last_col = function(x,d){
    x = cbind(x,matrix(x[,ncol(x)],nrow(x),ncol=d))
}

calc_rem_burden = function(prev_dt,agegroups_model,dE,dIp,dIs,dIa,dHosp,dDeath,ihr,ifr_dt){
    max_dates = prev_dt[,.(date=max(date)),by=.(country)]
    
    curr_prev_dt = merge(max_dates,prev_dt,by=c("country","date"))
    # Calculate rough estimate of number of days d until all susceptibles are
    # either infected or vaccinated if infection and vaccination rates remain 
    # constant at current levels
    curr_prev_dt[,d:=S/(nS_E+nS_V)]
    rem_days = curr_prev_dt[,.(d=max(d)),by=.(country)]
    
    countries = prev_dt[,unique(country)]
    ncountries = length(countries)
    
    proj_prev_list = vector("list",ncountries)
    rem_burden_list = vector("list",ncountries)
    na = length(agegroups_model)
    for (i in 1:ncountries){
        cntry = countries[i]
        d = rem_days[country==cntry,ceiling(d)]
        exposures = as.matrix(dcast(prev_dt[country==cntry],age_group_model ~ date,value.var = "exposures")[,-1])
        nt = ncol(exposures)
        # Repeat last column of exposures for d days, i.e. treat infection 
        # incidence as remaining constant
        exposures = rep_last_col(exposures,d)
        population = unique(prev_dt[country==cntry,.(age_group_model,population)])
        nS_V = as.matrix(dcast(prev_dt[country==cntry,.(date,age_group_model,nS_V)],age_group_model ~ date,value.var = "nS_V")[,-1])
        # Repeat last column of nS_V for d days, i.e. treat vaccination rate 
        # as remaining constant
        nS_V = rep_last_col(nS_V,d)
        q = as.matrix(dcast(prev_dt[country==cntry,.(date,age_group_model,q)],age_group_model ~ date,value.var = "q")[,-1])
        q = rep_last_col(q,d)
        r = as.matrix(dcast(prev_dt[country==cntry,.(date,age_group_model,r)],age_group_model ~ date,value.var = "r")[,-1])
        r = rep_last_col(r,d)
        s = as.matrix(dcast(prev_dt[country==cntry,.(date,age_group_model,s)],age_group_model ~ date,value.var = "s")[,-1])
        s = rep_last_col(s,d)
        
        Y = solve_diff_eqns(exposures,population[,population],nt+d,q,r,s,nS_V,covy,dE,dIp,dIs,dIa)
        
        indcs = reshape2::melt(Y[[1]])[,c("Var1","Var2")]
        age_grps = agegroups_model[indcs$Var1]
        dts = prev_dt[country==cntry,unique(date)][indcs$Var2[1:(na*nt)]]
        dts = c(dts,rep(seq.Date(dts[length(dts)]+1,dts[length(dts)]+d,by=1),each=na))
        for (j in seq_along(Y)){
            Y[[j]] = reshape2::melt(Y[[j]],value.name = names(Y)[j])[names(Y)[j]]
        }
        Y = do.call(cbind,Y)
        setDT(Y)
        Y[,`:=`(country=cntry,age_group_model=age_grps,date=dts)]
        
        # Merge with IHR and IFR data tables
        Y = merge(Y,ihr,by="age_group_model")
        # N.B. Change IFR so that it's the same for all countries with covidm 
        # age groups for consistency with covidm
        Y = merge(Y,ifr_dt,by=c("country","age_group_model"))
        # Merge time-varying vaccine efficacies against hospitalisation and death
        Y = merge(Y,prev_dt[country==cntry,.(country,age_group_model,date,eh)],by=c("country","age_group_model","date"),all.x=T)
        Y = merge(Y,prev_dt[country==cntry,.(country,age_group_model,date,em)],by=c("country","age_group_model","date"),all.x=T)
        # Forward fill NAs in vaccine efficacies, i.e. assume vaccine type 
        # proportions and variant proportions remain constant from now
        Y[,`:=`(eh=nafill(eh,"locf"),em=nafill(em,"locf")),by=.(age_group_model)]
        
        # Multiply infections by infection-hospitalisation rate and convolve 
        # infections with infection-to-hospitalisation delay to get 
        # hospitalisation curve and likewise for deaths with IFR
        Y[,hosp := ihr*(disc_conv(nS_E,dHosp)+(1-eh)*disc_conv(nV_E,dHosp)),by=.(age_group_model)]
        Y[,deaths := ifr*(disc_conv(nS_E,dDeath)+(1-em)*disc_conv(nV_E,dDeath)),by=.(age_group_model)]
        
        # Merge population
        Y = merge(Y,population,by="age_group_model")
        
        proj_prev_list[[i]] = Y
        
        # Sum hospitalisations and deaths from now until there are no more susceptibles
        rem_burden_list[[i]] = Y[date>=max(date)-d+1 & date<=max(date),.(population=mean(population),cum_hosp=sum(hosp),cum_deaths=sum(deaths)),by=.(age_group_model)]
        rem_burden_list[[i]][,`:=`(country=cntry,d=d)]
        cols = c("cum_hosp","cum_deaths")
        rem_burden_list[[i]][,(sub("cum","cum_inc",cols)):=lapply(.SD,function(x) x/population),.SDcols=cols]
    }
    proj_prev_dt = rbindlist(proj_prev_list)
    rem_burden_dt = rbindlist(rem_burden_list)
    rem_burden_dt[,age_group_model:=factor(age_group_model,levels = agegroups_model)]
    
    return(list(proj_prev_dt=proj_prev_dt,rem_burden_dt=rem_burden_dt))
}
