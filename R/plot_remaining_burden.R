plot_remaining_burden = function(rem_burden_dt,ovrl_rem_burden_dt,dir_fig="./figs/"){
    # Plot vaccine coverage pyramids
    p1 = ggplot(melt(rem_burden_dt,measure.vars = c("pop_prop_u","pop_prop_i","pop_prop_v")),aes(x=value,y=age_group_model,alpha=variable)) +
        geom_col(position = "stack") +
        scale_alpha_manual(name="",values=c(0.3,0.6,1),c("Unvaccinated &\nunexposed","Previously\ninfected","Partially/fully\nvaccinated &\nuninfected")) +
        labs(x="Proportion of population",y="Age group") +
        theme(axis.text.x=element_text(size=11)) +
        facet_wrap(~country)
    # ggsave(paste0(dir_fig,"vax_cov_pop_pyramids.png"),width = 10,height = 8)
    
    # Plot maximum remaining hospitalisations and deaths by age
    p2 = ggplot(melt(rem_burden_dt,measure.vars = c("cum_inc_hosp_u","cum_inc_hosp_i","cum_inc_hosp_v")),aes(y=age_group_model)) +
        geom_col(aes(x=value*1e5,alpha=variable),position="stack") +
        geom_errorbarh(aes(xmin=cum_inc_hosp_q95l*1e5,xmax=cum_inc_hosp_q95u*1e5),height=0.5) +
        scale_alpha_manual(name="",values=c(0.3,0.6,1),labels=c("Unvaccinated &\nunexposed","Previously\ninfected","Partially/fully\nvaccinated &\nuninfected")) +
        coord_cartesian(xlim = c(0,rem_burden_dt[!(country %in% c("Germany","Greece","Romania","Slovakia")),max(cum_inc_hosp_q95u)*1e5])) +
        labs(x="Maximum remaining hospitalisations/100,000 population",y="Age group") +
        facet_wrap(~country)
    # ggsave(paste0(dir_fig,"rem_hosps_by_age.png"),width = 10,height = 8)
    
    p3 = ggplot(melt(rem_burden_dt,measure.vars = c("cum_inc_deaths_u","cum_inc_deaths_i","cum_inc_deaths_v")),aes(y=age_group_model)) +
        geom_col(aes(x=value*1e5,alpha=variable),position="stack") +
        geom_errorbarh(aes(xmin=cum_inc_deaths_q95l*1e5,xmax=cum_inc_deaths_q95u*1e5),height=0.5) +
        scale_alpha_manual(name="",values=c(0.3,0.6,1),labels=c("Unvaccinated &\nunexposed","Previously\ninfected","Partially/fully\nvaccinated & uninfected")) +
        coord_cartesian(xlim = c(0,rem_burden_dt[!(country %in% c("Germany","Greece","Romania","Slovakia")),max(cum_inc_deaths_q95u)*1e5])) +
        labs(x="Maximum remaining deaths/100,000 population",y="Age group") +
        facet_wrap(~country) +
        theme(legend.position = "bottom")
    # ggsave(paste0(dir_fig,"rem_deaths_by_age.png"),width = 10,height = 8)
    
    # Plot maximum overall remaining burden:
    # hospitalisations against vaccine coverage
    p4 = ggplot(ovrl_rem_burden_dt,aes(x=cum_prop_v)) +
        geom_point(aes(y=cum_inc_hosp*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
        geom_errorbar(aes(ymin=cum_inc_hosp_q95l*1e5,ymax=cum_inc_hosp_q95u*1e5),width = 0.01,alpha=0.5) +
        geom_text(aes(y=cum_inc_hosp*1e5,label=country),size=3.5,hjust=-0.2,vjust=0.4) +
        xlim(NA,1) +
        labs(x="Proportion who have received at least one dose",
             y="Maximum remaining hospitalisations/100,000 population",
             size="Cumulative\nproportion\ninfected") +
        scale_y_log10()
    # ggsave(paste0(dir_fig,"rem_hosps_vs_prop_vax.png"),width = 8,height = 6.4)
    
    p = (p1 + theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1)) + p2 + theme(axis.text.x=element_text(angle=45,hjust=1))) / 
        p4 + plot_layout(heights = c(1,1.5)) + plot_annotation(tag_levels = 'A')
    ggsave(paste0(dir_fig,"vax_cov_and_rem_hosps.png"),plot=p,width = 16,height = 15)
    ggsave(paste0(dir_fig,"vax_cov_and_rem_hosps.pdf"),plot=p,width = 16,height = 15)
    
    # deaths against vaccine coverage
    p5 = ggplot(ovrl_rem_burden_dt,aes(x=cum_prop_v)) +
        geom_point(aes(y=cum_inc_deaths*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
        geom_errorbar(aes(ymin=cum_inc_deaths_q95l*1e5,ymax=cum_inc_deaths_q95u*1e5),width = 0.01,alpha=0.5) +
        geom_text(aes(y=cum_inc_deaths*1e5,label=country),size=4,hjust=-0.15,vjust=0.4) +
        xlim(NA,1) +
        labs(x="Proportion who have received at least one dose",
             y="Maximum remaining deaths/100,000 population",
             size="Cumulative proportion infected") +
        scale_y_log10() +
        theme(legend.position = "bottom")
    # ggsave(paste0(dir_fig,"rem_deaths_vs_prop_vax.png"),width = 8,height = 6.4)
    ggsave(paste0(dir_fig,"vax_cov_and_rem_deaths.png"),plot_grid(p3,p5,rel_widths=c(1,1),labels=c("A","B")),width = 14.4,height = 6)
    
    # hospitalisations against proportion aged 60+
    p6 = ggplot(ovrl_rem_burden_dt,aes(x=prop_pop_60plus)) +
        geom_point(aes(y=cum_inc_hosp*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
        geom_errorbar(aes(ymin=cum_inc_hosp_q95l*1e5,ymax=cum_inc_hosp_q95u*1e5),width = 0.002,alpha=0.5) +
        geom_text(aes(y=cum_inc_hosp*1e5,label=country),size=4,hjust=-0.15,vjust=0.4) +
        xlim(0.23,0.31) +
        labs(x="Proportion aged 60+ years",
             y="Maximum remaining hospitalisations/100,000 population",
             size="Cumulative proportion infected") +
        scale_y_log10()
    # ggsave(paste0(dir_fig,"rem_hosps_vs_prop_60plus.png"),width = 8,height = 6.4)
    
    # deaths against proportion aged 60+
    p7 = ggplot(ovrl_rem_burden_dt,aes(x=prop_pop_60plus)) +
        geom_point(aes(y=cum_inc_deaths*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
        geom_errorbar(aes(ymin=cum_inc_deaths_q95l*1e5,ymax=cum_inc_deaths_q95u*1e5),width = 0.002,alpha=0.5) +
        geom_text(aes(y=cum_inc_deaths*1e5,label=country),size=4,hjust=-0.15,vjust=0.4) +
        xlim(0.23,0.31) +
        labs(x="Proportion aged 60+ years",
             y="Maximum remaining deaths/100,000 population",
             size="Cumulative\nproportion\ninfected") +
        scale_y_log10()
    # ggsave(paste0(dir_fig,"rem_deaths_vs_prop_60plus.png"),width = 8,height = 6.4)
    
    p8 = plot_grid(p6+theme(legend.position="none"),
                   p7+theme(legend.position="none"),
                   labels=c("A","B"))
    l = get_legend(p6 + theme(legend.position = "bottom")) # + theme(legend.box.margin = margin(0,0,0,12)))
    ggsave(paste0(dir_fig,"rem_hosps_and_deaths_vs_prop_60plus.png"),plot_grid(p8,l,nrow=2,rel_heights = c(1,0.1)),width = 15,height = 6,bg = "white")
}
