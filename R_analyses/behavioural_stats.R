# Initialise ------

pacman::p_load('tidyverse','magrittr','patchwork','multcomp','emmeans','afex','ggpubr',
               'R.matlab')

# set random seed for reproducibility
set.seed(17)

# define custom functions

My_Theme = theme(
  axis.title.x = element_text(size = 25,vjust = -3),
  axis.text.x = element_text(size = 20,vjust = 0),#vjust adjust distance from tick to label
  axis.text.y = element_text(size = 20, margin = margin(t = 0, r = 5, b = 0, l = 0)),
  axis.title.y = element_text(size = 25,vjust = 3),
  axis.ticks = element_line(size=1),
  legend.text=element_text(size=25),
  legend.title = element_text(size=0),
  axis.ticks.length = unit(0.2, "cm"),
  plot.margin = margin(1, 1, 1,1, "cm"))


prep_long <- function(ctr,PD,dep.lab){

  ctr$group <- 'Ctr'
  ctr$ID <-seq(1,nrow(ctr))
  PD$group <-'PD'
  PD$ID <-seq(nrow(ctr)+1,nrow(ctr)+nrow(PD))
  allRT <-rbind(PD,ctr)
  names(allRT)<-c('Low_Low','High_Low','Low_High','High_High','Group','ID' )
  allRT_long <- tidyr::gather(allRT, Condition, RT, 'Low_Low':'High_High', factor_key=TRUE)
  allRT_long <- cbind(allRT_long,stringr::str_split_fixed(allRT_long$Condition,'_',2)) #split condition string column into multiple columns
  names(allRT_long)[5:6] <- c('PU','AU')
  plyr::rename(allRT_long,c('RT' = dep.lab))
  return(allRT_long)
}

# load data
m <- readMat(here::here('data','RTdataPD.mat'))
PD <- as.data.frame(m$ResPD[,,1])

# RT -----

#Mixed-design rm anova with group (ctr PD) as between factor and uncertainty (Perceptual Action) as within factor

ctr <- as.data.frame(m$ResCtr[,,1])

allRT_long <- prep_long(ctr,PD[,1:4],'RT')
allRT_long$PU %<>% factor(levels = c('Low','High'))
allRT_long$AU %<>% factor(levels = c('Low','High'))
allRT_long$ID %<>% factor()

# We will use logRT for stat but plot natural RT for readability
allRT_long$logRT <- log10(allRT_long$RT)
aov_ez('ID','logRT',allRT_long,within = c('PU','AU'), between = 'Group',factorize = FALSE)


allRT_long$AU <- factor(allRT_long$AU, levels = c('Low','High'))
allRT_long$PU <- factor(allRT_long$PU, levels = c('Low','High'))
allRT_long$Group%<>% factor()

AU.labs <- c('Low Action \nUncertainty','High Action \nUncertainty');
names(AU.labs) <- c('Low','High')

# Plot RTs - Figure 3, Panel B ----
ggplot(allRT_long, aes(x=PU, y=RT, fill=Group)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  facet_wrap(~AU,labeller=labeller(AU = AU.labs)) +
  scale_fill_manual(values = c('#B2182B','#2166AC'))+
  labs( x = "Perceptual Uncertainty",y="RT (s)")+theme_minimal()+My_Theme+theme(strip.text = element_text(size = 25))
  


#Bayesian stats
bf_full<-BayesFactor::generalTestBF(RT~AU*PU*Group+ID, whichRandom = 'ID',data=allRT_long,neverExclude = 'ID')
bf_null<-BayesFactor::generalTestBF(RT~ID, whichRandom = 'ID',data=allRT_long)
bf_final<- bf_full/bf_null
bayestestR::bayesfactor_inclusion(bf_final,match_models=TRUE)

#Plot densities for figure 1B

dfplot<-allRT_long 
 levels(dfplot$PU)<-c('Low PU', 'High PU')
 levels(dfplot$AU)<-c('Low AU', 'High AU')
 
  ggplot(dfplot,aes(x=RT,fill = Group)) + geom_density(alpha = 0.8,show.legend = TRUE)+
   facet_grid(PU~AU)+scale_x_continuous(name = 'Reaction Time (secs)',limits = c(0.2,1.7),expand = c(0,0))+
   scale_y_continuous(name = "Density", expand = c(0, 0))+
   #scale_fill_manual(values = c('#B2182B','#2166AC'))+
    scale_fill_manual(values = c('#E41A1C',"#377EB8"))+
  theme_minimal()+theme(
   #legend.position = "none",
   #plot.title = element_text(size = 16, face = "bold", family = "Helvetica", colour = "steelblue4", vjust = 4, hjust = 0.5),
   #plot.subtitle = element_text(size = 14,  family = "Helvetica",colour = "steelblue4", vjust = 2, hjust = 0.5),
   axis.title.x = element_text(size = 14, face = "bold", family = "Helvetica", vjust = -2),
   axis.title.y = element_text(size = 14, face = "bold", family = "Helvetica",  vjust = 2),
   axis.text.x = element_text(size = 12, family = "Helvetica"),
   axis.text.y =element_text(size = 12, family = "Helvetica"),
   panel.grid.minor = element_blank(),
   panel.grid.major = element_line(size = 0.4),
     strip.text = element_text(size = 12, family = "Helvetica"),
   plot.margin = unit(c(1,1,1,1), "cm"),
   legend.position = "bottom",
   legend.justification = "right",
   legend.margin = margin(15, 0, 1.5, 0, "pt"),
   legend.spacing.x = grid::unit(3, "pt"),
   legend.spacing.y = grid::unit(0, "pt"),
   legend.box.spacing = grid::unit(0, "pt"),
   legend.title=element_blank(),
   legend.text=element_text(size = 12, family = "Helvetica")
 )

  
 



# Accuracy --------
Ctracc <- as.data.frame(m$ResCtr[,,2])

PDacc <- as.data.frame(m$ResPD[,,2])
Ctracc$ID <- seq(1,21)
Ctracc$Group <- 'Ctr'
PDacc$Group <- 'PD'
PDacc$ID <- seq(22,38)

acc <- rbind(Ctracc,PDacc)

#afex requires long format
names(acc)[1:4] <- c('Low_Low','High_Low','Low_High','High_High')
acc_long <- tidyr::gather(acc,Condition,Accuracy,'Low_Low':'High_High')
acc_long <- cbind(acc_long,stringr::str_split_fixed(acc_long$Condition,'_',2)) #split condition string column into multiple columns
names(acc_long)[5:6] <- c('PU','AU')
acc_long$PU <- factor(acc_long$PU,levels = c('Low','High'))
acc_long$AU <- factor(acc_long$AU,levels = c('Low','High'))
acc_long$ID %<>% factor()

PD$UPDRS <-scale(t(m$UPDRS))
PD$LEDD  <-scale(t(m$LEDD))

#Mixed-design rm anova with group (ctr PD) as between factor and uncertainty (Perceptual Action) as within factor
(acc_mod <- aov_ez('ID','Accuracy',acc_long,within = c('PU','AU'),between = 'Group'))

#multiplicity adjusted p values for each linear hypotheses controlling the family-wise error rate (FWER)
#free: multiple testing procedures under free combinations (Westfall et al., 1999)
#PD Committed more errors
main_acc <- emmeans(acc_mod,~Group)
summary(as.glht(pairs(main_acc)),test=adjusted('bonferroni'))
#more errors for High PU
main_acc <- emmeans(acc_mod,~PU)
summary(as.glht(pairs(main_acc)),test=adjusted('bonferroni'))
#more errors for Low AU
main_acc <- emmeans(acc_mod,~AU)
summary(as.glht(pairs(main_acc)),test=adjusted('bonferroni'))

(inter_acc <- emmeans(acc_mod,~AU|Group))
pairs(inter_acc)
summary(as.glht(pairs(inter_acc)),test=adjusted('bonferroni'))

# Prepare for plots

acc_long$AU <-factor(acc_long$AU, levels = c('Low','High'))
AU.labs <- c('Low Action \nUncertainty','High Action \nUncertainty');
names(AU.labs) <- c('Low','High')


# Plot accuracy - Figure 3, Panel C ----

ggplot(acc_long, aes(x=PU, y=Accuracy, fill=Group)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  facet_wrap(~AU,labeller=labeller(AU = AU.labs)) +
  scale_fill_manual(values = c('#B2182B','#2166AC'))+
  labs( x = "Perceptual Uncertainty",y="Accuracy p(c)/chance")+theme_minimal()+My_Theme+theme(strip.text = element_text(size = 25))


ggplot(acc_long, aes(x=Group, y=Accuracy, fill=Group)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_fill_manual(values = c('#B2182B','#2166AC'))+
  labs( x = "Group",y="Accuracy p(c)/chance")+theme_minimal()+My_Theme+theme(strip.text = element_text(size = 25))


acc_long$ID %<>% factor()
(acc_mod <- mixed(Accuracy~AU*PU*Group+(1|ID),acc_long))

(inter_acc <- emmeans(acc_mod,~AU|Group))
summary(as.glht(pairs(inter_acc)),test=adjusted('bonferroni'))
#Bayesian stats
bf_full<-BayesFactor::generalTestBF(Accuracy~AU*PU*Group+ID, whichRandom = 'ID',data=acc_long,neverExclude = 'ID')
bf_null<-BayesFactor::generalTestBF(Accuracy~ID, whichRandom = 'ID',data=acc_long)
bf_final<- bf_full/bf_null
bayestestR::bayesfactor_inclusion(bf_final,match_models=TRUE)

# Thresholds -----

th<- as.data.frame(readMat(here::here('data','Thr_group.mat')))
th$note <- NULL
th$Group <- ifelse(th$groupTh.2==1,'Ctr','PD')
th$groupTh.2 <- NULL
names(th)[1:5]<- c('ID','low_low','high_low','low_high','high_high')

#long format for rm-ANOVA and plotting 
long_th <- tidyr::gather(th,condition,th, 'low_low':'high_high')
long_th <- cbind(long_th,stringr::str_split_fixed(long_th$condition,'_',2)) #split condition string column into multiple columns
names(long_th)[5:6] <- c('PU','AU')
long_th$AU    %<>% factor( levels = c('low','high'))
long_th$PU    %<>% factor( levels = c('low','high'))
long_th$ID    %<>% factor()
long_th$Group %<>% factor()


#Mixed-design rm anova with group (ctr PD) as between factor and uncertainty (Perceptual Action) as within factor
aov_ez('ID','th',long_th,within = c('PU','AU'), between = 'Group',factorize = FALSE)

#Bayesian stats
bf_full<-BayesFactor::generalTestBF(th~AU*PU*Group+ID, whichRandom = 'ID',data=long_th,neverExclude = 'ID')
bf_null<-BayesFactor::generalTestBF(th~ID, whichRandom = 'ID',data=long_th)
bf_final<- bf_full/bf_null
bayestestR::bayesfactor_inclusion(bf_final,match_models=TRUE)


# Prepare for plots
AU.labs <- c('Low Action \nUncertainty','High Action \nUncertainty');
names(AU.labs) <- c('low','high')


# Plot thresholds - Figure 3, Panel A ----

ggplot(long_th, aes(x=PU, y=th, fill=Group)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  facet_wrap(~AU,labeller=labeller(AU = AU.labs)) +
  scale_fill_manual(values = c('#B2182B','#2166AC'))+
  labs( x = "Perceptual Uncertainty",y="Coherence (prop)")+theme_minimal()+My_Theme+theme(strip.text = element_text(size = 25))

