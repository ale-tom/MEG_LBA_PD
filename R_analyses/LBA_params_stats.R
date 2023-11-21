library(R.matlab)
library(tidyverse)
library(magrittr)
library(emmeans)
library(patchwork)
library(multcomp)
library(bbplot)

set.seed(17)

# ---- Load and prepare the data ----

m <- readMat(here::here('data','PD_model_pars.mat'))
m = m$cellV
m = rbind(m[,,1],m[,,2])
m <- m[complete.cases(m),]#remove NaNs

df_lba_par = as.data.frame(m)
df_lba_par = data.frame(col = unlist(df_lba_par,use.names = FALSE))
samplesize = nrow(m)
nctr = 21; npd = 17;
names(df_lba_par)[1] <- 'V'
df_lba_par$PID <- as.factor(rep(1:samplesize,4))
df_lba_par$PU  <- as.factor(rep(c(rep('low',samplesize),rep('high',samplesize)),2))
df_lba_par$AU  <- as.factor(c(rep('low',samplesize*2),rep('high',samplesize*2)))
df_lba_par$Group  <- as.factor(c(rep(c(rep('Ctr',nctr),rep('PD',npd)),4)))


#reorder columns
df_lba_par <- df_lba_par[,c('PID','PU','AU','Group','V')]
#sort data
df_lba_par <- df_lba_par[order(df_lba_par$PID), ]
df_lba_par <- within(df_lba_par, {
  PID   <- factor(PID)
  PU <- factor(PU,levels = c('low','high'))
  AU <- factor(AU,levels = c('low','high'))
  Group  <- factor(Group)
})

# ----Mixed-design repeted-measures ANOVA on LBA parameters----

# Mixed-design repeated measures ANOVA - Frequentist & Bayesian
(v_mod<-afex::aov_ez('PID','V',df_lba_par,within = c('PU','AU'), between = 'Group',factorize = TRUE))
bf_v <- BayesFactor::anovaBF(V~Group*AU*PU + PID,whichRandom = 'PID',data = df_lba_par,)

##Bayes interaction AU Group (note, this will slightly vary each time)
bf_au_grp <- bf_v[17]/bf_v[12]


# Post-hoc tests for AU*Group interaction
(inter_AU_group <- emmeans(v_mod,~AU*Group))
int_au<-summary(as.glht(pairs(inter_AU_group)),test=adjusted('bonferroni'))


# Define the contrasts
contrasts_list <- list(
  interaction = c(1, -1, -1, 1)
)

# Apply the contrasts to the emmeans
interaction_contrast_result <- contrast(inter_AU_group, method = contrasts_list)

# Test the interaction contrast with adjustment for multiple comparisons
summary(interaction_contrast_result, test=adjusted('bonferroni'))

#contrast    estimate    SE df t.ratio p.value
#interaction     1.18 0.498 36 2.377   0.0229 



# ---- Figure 7 (right panel) ------

# Prepare for plots

(inter_PU_group <- emmeans(v_mod,~PU*Group))
int_pu<-summary(as.glht(pairs(inter_PU_group)),test=adjusted('bonferroni'))



df_bar_au<-data.frame(matrix(nrow=2))
df_bar_au$Est<-c(int_au$test$coefficients[1][[1]], int_au$test$coefficients[6][[1]])
df_bar_au$se <-c(int_au$test$sigma[1][[1]], int_au$test$sigma[6][[1]])
df_bar_au$Group<-c('Ctr','PD')

df_bar_pu<-data.frame(matrix(nrow=2))
df_bar_pu$Est<-c(int_pu$test$coefficients[1][[1]], int_pu$test$coefficients[6][[1]])
df_bar_pu$se <-c(int_pu$test$sigma[1][[1]], int_pu$test$sigma[6][[1]])
df_bar_pu$Group<-c('Ctr','PD')

#Results are averaged over the levels of: PU and 
df_bar_pu<-data.frame(matrix(nrow=2))
df_bar_pu$Est<-c(int_pu$test$coefficients[1][[1]], int_pu$test$coefficients[6][[1]])
df_bar_pu$se <-c(int_pu$test$sigma[1][[1]], int_pu$test$sigma[6][[1]])
df_bar_pu$Group<-c('Ctr','PD')


scaleFUN <- function(x) sprintf("%.1f", x)
create_plot_data_v2 <- function(data, df_bar_au, grouping_var) {
  
  # Ensure the grouping variable is either 'AU' or 'PU'
  if (!grouping_var %in% c("AU", "PU")) {
    stop("grouping_var must be either 'AU' or 'PU'")
  }
  
  # Compute the individual differences
  difference_data <- data %>%
    group_by(Group, PID) %>%
    summarise(diff = mean(V[!!sym(grouping_var) == "low"]) - mean(V[!!sym(grouping_var) == "high"])) %>%
    ungroup()
  
  # Create a dataframe with unique PIDs and their respective groups
  pid_data <- data.frame(PID = unique(data$PID), Group = rep(unique(data$Group), times = c(sum(data$Group == "Ctr"), sum(data$Group == "PD"))))
  
  # Join the df_bar_au with the pid_data
  expanded_dfbar <- left_join(pid_data, df_bar_au, by = "Group")
  
  # Bind the difference data columns to the expanded df_bar_au
  plot_data <- left_join(expanded_dfbar, difference_data, by = c("Group", "PID"))
  
  return(plot_data %>% distinct())
}


plot_data_AU = create_plot_data_v2(df_lba_par, df_bar_au, 'AU')
plot_data_PU = create_plot_data_v2(df_lba_par, df_bar_pu, 'PU')

plot_AU<-ggplot(data=plot_data_AU, aes(x=Group, y=Est,fill=Group)) +
  geom_bar(stat="identity", color='Black', position=position_dodge(),show.legend = FALSE)+
  geom_point(aes(y = diff), position =position_jitterdodge(),alpha=0.3, show.legend = FALSE) +  # Add individual points
  coord_cartesian(ylim = c(0, 8))+
  labs(title="",subtitle = 'Low-High AU', x="Group", y = "Low - High AU")+
  theme_minimal()+bbc_style()+
  scale_fill_manual(values=c('#B2182B','#2166AC'))+scale_y_continuous(labels=scaleFUN)

plot_PU<-ggplot(data=plot_data_PU, aes(x=Group, y=Est,fill=Group)) +
  geom_bar(stat="identity", color='Black', position=position_dodge(),show.legend = FALSE)+
  geom_point(aes(y = diff), position =position_jitterdodge(),alpha=0.3, show.legend = FALSE) +  # Add individual points
  coord_cartesian(ylim = c(0, 8))+
  labs(title='Accumulation rate',subtitle = 'Low-High PU', x="Group", y = "Low - High PU")+
  theme_minimal()+bbc_style()+
  scale_fill_manual(values=c('#B2182B','#2166AC'))+scale_y_continuous(labels=scaleFUN)

plot_PU + plot_AU


sjPlot::save_plot(here::here('Model','figures','AccrateInteractions_single_points.jpg'),fig=last_plot(),width =13, height = 10, dpi = 300)
dev.off()
