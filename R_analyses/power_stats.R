
pacman::p_load('R.matlab', 'magrittr','lme4','lmerTest','emmeans','tidyverse','bbplot','patchwork')

set.seed(17)

# Load and prepare the data ----

# m is 2 x 3D structure (Nsubj x 96rois x 4 conditions)
m <- readMat(here::here('data','beta_meanpow_roi_stat.mat'))

roi_table <- as.data.frame(read.csv(here::here('data','ROI_table.csv')))
roi_table$ROI <- factor(roi_table$ROI)
levels(roi_table$ROI)[18]<-'Supplementary Motor Cortex'
levels(roi_table$ROI)[10]<- 'Heschls Gyrus'
roi_table %<>% rbind(roi_table)  
roi_table$ROInum <- 1:96
names(roi_table)[[5]]<-'ROI_Labels'



# transform into a 2D array, by concatenating the 3rd dimension along the 1st (i.e. row-wise) + add column corresponding to condition and subject id
ctr <- m$H.cell[[1]][[1]]

ctr<-Reduce(rbind, sapply(1:dim(ctr)[3], function(i) {
  x <- data.frame(ctr[,,i])
  x$condition <- i
  x$subject <- seq(1,nrow(ctr[,,1]))
  x 
  }, simplify = FALSE))

ctr$Group <- 'Ctr'


pd <- m$H.cell[[2]][[1]]

pd<-Reduce(rbind, sapply(1:dim(pd)[3], function(i) {
  x <- data.frame(pd[,,i])
  x$condition <- i
  x$subject <- seq(22,21+nrow(pd[,,1]))
  x 
}, simplify = FALSE))

pd$Group <- 'PD'


# concatenate ctr and pd
all_pow <- rbind(ctr,pd)

# name levels
all_pow$condition %<>% plyr::mapvalues(
                     from = seq(1,4),
                     to = c("low_low", "high_low", "low_high", "high_high"))

# rename ROIs recursively
all_pow %<>% dplyr::rename_at(1:96,~paste0("ROI_", seq(1,96)))

dorsal_path <- c(1, 3, 4, 7, 26, 17, 18, 31, 32, 23, 22, 24, 48)
dorsal_path <- c(dorsal_path, dorsal_path + 48)


# Moderation effect of PD on beta desynchronization ----

#lmer requires long format
long_pow <- tidyr::gather(all_pow,ROI,Amplitude, 'ROI_1':'ROI_96')
long_pow <- cbind(long_pow,stringr::str_split_fixed(long_pow$condition,'_',2)) #split condition string column into multiple columns
names(long_pow)[6:7] <- c('PU','AU')
long_pow$condition <- NULL
long_pow$subject<-as.factor(long_pow$subject)
long_pow$Group %<>% factor(levels = c('Ctr','PD'))
long_pow$AU %<>% factor(levels = c('low','high'))
long_pow$AU %<>% relevel(ref="low")
long_pow$PU %<>% factor(levels = c('low','high'))
long_pow$AU %<>% relevel(ref="low")

long_pow$ROI %<>% factor()

# Linear mixed effect model
beta_pow<-lmer(Amplitude~PU+AU+Group+PU*Group+AU*Group+(1|subject/ROI),data = long_pow,REML = F)

report::report(beta_pow)

# Post-hoc test
em_A_Beta <- emmeans::emmeans(beta_pow, ~Group*AU,lmerTest.limit=14592)
pairs(emmeans(em_A_Beta, ~ AU | Group), adjust='bonferroni')



# Reactivity: Interaction contrast between low AU-high AU differences between groups ----
# Define the contrasts
contrasts_list <- list(
  interaction = c(1, -1, -1, 1)
)

# Apply the contrasts to the emmeans
interaction_contrast_result <- contrast(em_A_Beta, method = contrasts_list)

# Test the interaction contrast with adjustment for multiple comparisons
summary(interaction_contrast_result, test=adjusted('bonferroni'))
#contrast    estimate    SE  df z.ratio p.value
#interaction    0.425 0.124 10944 3.429   0.0006 


# Plot reactivity interactions -----

create_plot_data <- function(data, dfbar, model_fit, grouping_var) {
  
  # Ensure the grouping variable is either 'AU' or 'PU'
  if (!grouping_var %in% c("AU", "PU")) {
    stop("grouping_var must be either 'AU' or 'PU'")
  }
  # Create a new data frame with every combination of `subject`, `AU`, `Group`, `PU`, and `ROI`
  new_data <- expand.grid(subject = unique(data$subject),
                          AU = unique(data$AU),
                          Group = ifelse(as.numeric(as.character(unique(data$subject))) <= 21, "Ctr", "PD"),
                          PU = unique(data$PU),  
                          ROI = unique(data$ROI))
  
  # Generate predictions
  new_data$predicted_values <- predict(model_fit, newdata = new_data)
  
  
  # Z-score scaling the predictions:
  data <- data %>%
          group_by(!!sym(grouping_var), Group) %>%
          mutate(z_amplitude = scale(Amplitude)) %>%
          ungroup()
  
  # Calculate the marginal mean for each level in the predictions:
  marginal_means <- new_data %>%
                    group_by(!!sym(grouping_var), Group) %>%
                    summarize(marginal_mean = mean(predicted_values)) %>%
                    ungroup()
  
  
  # Shift the Z-scored raw data to center it around the marginal mean for each level:
  data <- data %>%
          left_join(marginal_means, by = c(grouping_var, "Group")) %>%
          mutate(z_amplitude_shifted = z_amplitude + marginal_mean) %>%
          select(-marginal_mean)
  
  # Aggregate across ROIs
  mean_vals <- data %>%
               group_by(Group, subject, ROI) %>%
               summarise(diff = median(z_amplitude_shifted[!!sym(grouping_var) == "low"]) - median(z_amplitude_shifted[!!sym(grouping_var) == "high"])) %>%
               ungroup() %>% 
               group_by(Group,subject) %>%  
               summarise(avg_diff = median(diff))
  
  plot_data <- left_join(dfbar, mean_vals, by = "Group")
  
  return(plot_data)
}


# estimate coefficients interaction GROUP x Perceptual uncertainty
inter_au <- emmeans(beta_pow,~AU*Group)
int_au<-summary(as.glht(pairs(inter_au)),adjsted = 'bonferroni')
inter_pu <- emmeans(beta_pow,~PU*Group)
int_pu<-summary(as.glht(pairs(inter_pu)),adjust='bonferroni')



dfbar_au<-data.frame(matrix(nrow=2))
dfbar_au$Est<-c(int_au$test$coefficients[1][[1]], int_au$test$coefficients[6][[1]])
dfbar_au$se <-c(int_au$test$sigma[1][[1]], int_au$test$sigma[6][[1]])
dfbar_au$Group<-c('Ctr','PD')

dfbar_pu<-data.frame(matrix(nrow=2))
dfbar_pu$Est<-c(int_pu$test$coefficients[1][[1]], int_pu$test$coefficients[6][[1]])
dfbar_pu$se <-c(int_pu$test$sigma[1][[1]], int_pu$test$sigma[6][[1]])
dfbar_pu$Group<-c('Ctr','PD')


plot_data_AU<-create_plot_data(long_pow, dfbar_au, beta_pow, 'AU')
plot_data_PU<-create_plot_data(long_pow, dfbar_pu, beta_pow, 'PU')



# Plot for AU
plot_au <-ggplot(data = plot_data_AU, aes(x = Group, y = Est, fill = Group)) +
          geom_bar(stat = "identity", color = 'Black', 
                   position = position_dodge(0.9), show.legend = FALSE) +
          geom_point(aes(y = avg_diff), 
                     position =position_jitterdodge(),alpha=0.3, show.legend = FALSE) +  # Add individual points
          coord_cartesian(ylim = c(1, 3.5)) +
          labs(title = "", subtitle = 'Low-High AU', x = "Group", y = "Low - High AU") +
          theme_minimal() + 
          bbc_style() +
          scale_fill_manual(values = c('#B2182B', '#2166AC'))

# Plot for PU
plot_pu<-ggplot(data = plot_data_PU, aes(x = Group, y = Est, fill = Group)) +
         geom_bar(stat = "identity", color = 'Black', 
                  position = position_dodge(0.9), show.legend = FALSE) +
         geom_point(aes(y = avg_diff), 
                    position =position_jitterdodge(),alpha=0.3, show.legend = FALSE) +  # Add individual points
         coord_cartesian(ylim = c(-1, -3.5)) +
         labs(title='Power (Beta)',subtitle = 'Low-High PU', x="Group", y = "Low - High PU")+
         theme_minimal() + 
         bbc_style() +
         scale_fill_manual(values = c('#B2182B', '#2166AC'))

plot_pu + plot_au



# Differences in beta power between groups ----

# test whether beta desynchronization is diminished in PD vs controls
mean_pow <- long_pow %>% 
            mutate_at('ROI', str_remove,'ROI_') %>% filter(ROI %in% dorsal_path) %>%
            group_by(Group,subject)  %>% 
            summarize(MeanAmplitude = mean(Amplitude)) %>% 
            ungroup() 


# Wilcoxon rank sum test: beta desyncrhonization is significantly reduced in PD patients
mean_pow %>% select(Group,MeanAmplitude) %>% 
             gather(key = variable, value = value,-Group) %>%
             group_by(Group,variable) %>% 
             summarise(value = list(value)) %>%
             spread(Group,value) %>% 
             mutate(p_value = wilcox.test(unlist(Ctr),unlist(PD),paired = FALSE,exact = FALSE,alternative = 'less')$p.value,
                    W_value = wilcox.test(unlist(Ctr),unlist(PD),paired = FALSE,exact = FALSE,alternative = 'less')$statistic[[1]],
                    Ctr = mean(unlist(Ctr)),
                    PD = mean(unlist(PD)))

#variable        Ctr    PD p_value W_value
#<chr>         <dbl> <dbl>   <dbl>   <dbl>
# 1 MeanAmplitude -12.9 -9.75  0.0416     119

