#load packages
library(mice)
library(dplyr)

#set seeds and params
set.seed(4242)
repetitions = 2 #(number minus 1)
iter = 25 #iterations for MI

# generate complete data frame

data <- data.frame(X1=rnorm(1000), X2=rnorm(1000))
data$Y <- with(data, X1 + 4*X2 + rnorm(1000, sd=0.1))
write.csv(data, "data.csv", row.names=FALSE)

#extract summaries
a <- summary(lm(Y~X1+X2, data=data))$coefficients
colnames(a) = c("estimate", "std.error", "statistic", "p.value")
confo <- confint(lm(Y~X1+X2, data=data))
a <- cbind (a, confo)

#!!!!!repat from here based on the iter number (minus one)!!!!!#

#generate missingness with mar, mcar
data_mar10 = ampute(data, prop = 0.10, mech = "MAR")
data_mar25 = ampute(data, prop = 0.25, mech = "MAR")
data_mar40 = ampute(data, prop = 0.40, mech = "MAR")
data_mar55 = ampute(data, prop = 0.55, mech = "MAR")


data_mcar10 = ampute(data, prop = 0.10, mech = "MCAR")
data_mcar25 = ampute(data, prop = 0.25, mech = "MCAR")
data_mcar40 = ampute(data, prop = 0.40, mech = "MCAR")
data_mcar55 = ampute(data, prop = 0.55, mech = "MCAR")


#create empty dfs for summaries
Pattern_list <- list(mget(ls(pattern = "data_")))[[1]]
mod_summaries <- data.frame()
mod_summaries_sd <- list()

#complete case analyses
for (i in Pattern_list) {
  sumX <- as.data.frame(summary(lm(Y~X1+X2, data=i$amp))$coefficients) #get summary
  colnames(sumX) = c("estimate", "std.error", "statistic", "p.value")
  conf <- confint(lm(Y~X1+X2, data=i$amp))
  sumX <- cbind(sumX, conf)
  name1 = paste0(i$mech, i$prop, "cc")
  sumX <- cbind(sumX, name1)
  mod_summaries <- bind_rows(list(mod_summaries, sumX)) #combine summary + sds for output
}

#impute

#select imp method
method = "pmm"

#impute N times
x <- 1
repeat {
  print(x)

mimp_mar10 <- mice(data_mar10$amp, maxit = iter, m = 10, meth = method)
mimp_mar25<- mice(data_mar25$amp, maxit = iter, m = 25, meth = method)
mimp_mar40 <- mice(data_mar40$amp, maxit = iter, m = 40, meth = method)
mimp_mar55 <- mice(data_mar55$amp, maxit = iter, m = 55, meth = method)

mimp_mcar10 <- mice(data_mcar10$amp, maxit = iter, m = 10, meth = method)
mimp_mcar25 <- mice(data_mcar25$amp, maxit = iter, m = 25, meth = method)
mimp_mcar40 <- mice(data_mcar40$amp, maxit = iter, m = 40, meth = method)
mimp_mcar55 <- mice(data_mcar55$amp, maxit = iter, m = 55, meth = method)

#get the list of imputed sets, analayse and extract summaries
Pattern1_list <- list(mget(ls(pattern = "mimp")))[[1]]

for (i in Pattern1_list) {
  fit <- with(data=i,exp=lm(Y~X1+X2)) #linear model
  sumX <- summary(pool(fit), conf.int = T) #pool results
  name = as.data.frame(as.character(i$call$data)[2])
  name1 = paste0(name[,1],"mi")
  sumX <- bind_cols(sumX, as.data.frame(name1))#intialise name
  mod_summaries <- bind_rows(list(mod_summaries, sumX))
}

x = x+1
if (x == repetitions){
  break
}
}

#single impute and analyse

x <- 1
repeat {
  print(x)

simp_mar10 <- mice(data_mar10$amp, maxit = iter, m = 1, meth = method)
simp_mar25 <- mice(data_mar25$amp, maxit = iter, m = 1, meth = method)
simp_mar40 <- mice(data_mar40$amp, maxit = iter, m = 1, meth = method)
simp_mar55 <- mice(data_mar55$amp, maxit = iter, m = 1, meth = method)

simp_mcar10 <- mice(data_mcar10$amp, maxit = iter, m = 1, meth = method)
simp_mcar25 <- mice(data_mcar25$amp, maxit = iter, m = 1, meth = method)
simp_mcar40 <- mice(data_mcar40$amp, maxit = iter, m = 1, meth = method)
simp_mcar55 <- mice(data_mcar55$amp, maxit = iter, m = 1, meth = method)


Pattern2_list <- list(mget(ls(pattern = "simp")))[[1]]

  for (i in Pattern2_list) {
    sumX <- as.data.frame(summary(lm(Y~X1+X2, data=complete(i)))$coefficients) #get summary
    colnames(sumX) = c("estimate", "std.error", "statistic", "p.value")
    confs <- confint(lm(Y~X1+X2, data=complete(i)))
    sumX <- cbind(sumX, confs)
    name = as.data.frame(as.character(i$call$data))
    name1 = paste0(name[2,],"si")
    sumX <- cbind(sumX, name1)
    sumX$term = rownames(sumX)
    mod_summaries <- bind_rows(list(mod_summaries, sumX)) #combine summary + sds for output
  }

x = x+1
if (x == repetitions){
  break
}
}

#combine the cummaries and add the full set summary at the bottom
mod_summaries <- bind_rows(list(mod_summaries, as.data.frame(a)))
mod_summaries$name1[nrow(mod_summaries)] = c("original")
mod_summaries$name1[nrow(mod_summaries)-1] = c("original")
mod_summaries$name1[nrow(mod_summaries)-2] = c("original")
summaries <- mod_summaries #keep a backup in case


#calculate the bias
#check this after changing numer of predictors!!!!

#remove the intercep
toDelete <- seq(1, nrow(mod_summaries), 3)
mod_summaries = mod_summaries[ -toDelete,]

#create a sequence (even nrs X1, odd nrs X2)
biasX1 <- seq(1, nrow(mod_summaries), 2)
biasX2 <- seq(0, nrow(mod_summaries), 2)
biasX2 = biasX2[-1] #remove the 0

#create dfs to store X1 and X2 results separately
result1 = data.frame()
result2 = data.frame()
resultSE = data.frame()

#get se from the full model to calculate bias
result1se = mod_summaries[biasX1,2]
result2se = mod_summaries[biasX2,2]

#get the X1 and X2 estimate to subtract from the estimate from imputation
X1basetoremove = mod_summaries[c(nrow(mod_summaries)-1),c(1,3:6)]
X2basetoremove = mod_summaries[nrow(mod_summaries),c(1,3:6)]

#biasX1
for(i in biasX1){
  resultX1 = mod_summaries[i,c(1,3:6)] - X1basetoremove
  result1 = rbind(result1, resultX1)
}

#biasX2
for(i in biasX2){
  resultX2 = mod_summaries[i,c(1,3:6)] - X2basetoremove
  result2 = rbind(result2, resultX2)
}

#add the names
result1$names = mod_summaries$name1[biasX1]
result2$names = mod_summaries$name1[biasX2]

#bind in one df
result = rbind(result1, result2)

#calculate bias for SE (Pythagorean distance)
seX1 = mod_summaries[c(nrow(mod_summaries)-1),2]
seX2 = mod_summaries[nrow(mod_summaries),2]

result1se = mod_summaries[biasX1,2]
result2se = mod_summaries[biasX2,2]

for(i in 1: length(result1se)){
  resultSE1 = (sqrt(result1se[i]^2 + result2se[i]^2))/sqrt(seX1^2 + seX2^2)
  resultSE = rbind(resultSE, resultSE1)
}

#add the std error twice (for X1 and X2) - it's the same
result$std.error = rep(resultSE[,1], 2)
  

# PLOTTING
require(cowplot)
require(ggplot2)

#order cats
result$names <- factor(result$names, levels = c("original", "MAR0.1cc", "data_mar10mi", "data_mar10si", "MAR0.25cc", "data_mar25mi", "data_mar25si", "MAR0.4cc", "data_mar40mi", "data_mar40si", "MAR0.55cc", "data_mar55mi", "data_mar55si", "MCAR0.1cc", "data_mcar10mi", "data_mcar10si", "MCAR0.25cc", "data_mcar25mi", "data_mcar25si", "MCAR0.4cc", "data_mcar40mi", "data_mcar40si", "MCAR0.55cc", "data_mcar55mi", "data_mcar55si"))

#set theme for ggplot (will be overriden later)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)

#aggregate the data dropping the full dataset level and setting it as 0 or 1 for SE
res <- aggregate(cbind(`2.5 %`, `97.5 %`, estimate, std.error)~names,data=result,FUN=mean)
res = res[-1,] #remove the full datatset results
#add mechanism, number of iterations and proper names as separate columns
res$mech = rep(c("MAR", "MCAR"), times = c(12, 12))
res$num = rep(c("10", "25", "40", "55"), times = c(3, 3, 3, 3))
res$names = c("MAR CC", "MAR MI", "MAR SI", "MAR CC", "MAR MI", "MAR SI", "MAR CC", "MAR MI", "MAR SI", "MAR CC", "MAR MI", "MAR SI", "MCAR CC", "MCAR MI", "MCAR SI", "MCAR CC", "MCAR MI", "MCAR SI", "MCAR CC", "MCAR MI", "MCAR SI", "MCAR CC", "MCAR MI", "MCAR SI")
res$Type = rep(c("Complete case", "Multiple imputation", "Single imputation"), times = 8)

#set colours for plotting
cols = rep(c("black", "red2", "dodgerblue3"), times = 8)

#plot estimate
asd <- ggplot(res, aes(y = names, x = estimate)) +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                color = cols, alpha = 0.7, linewidth = 0.7) +
  geom_vline(aes(xintercept = 0), color = "red", 
             linewidth = 0.5, alpha = 0.8, linetype = "dashed") +  
  geom_point(color = cols, size = 2) +
  facet_grid(num ~ ., scales = "free", space = "free") +
  theme_light() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = cols),
        axis.title = element_text(size = 15,
                                  color = "black",
                                  face = "bold")) + 
  labs(x = "Estimate", y = "Missingness pattern") +
  coord_cartesian(xlim = c(-0.035, 0.035))

#plot SE
asd1 <- ggplot(res, aes(y = std.error, x = names)) +
  geom_segment( aes(x = names, xend=names, y=0, yend=std.error), color=cols) +
  geom_hline(aes(yintercept = 1), color = "red", 
             linewidth = 0.7, alpha = 0.8, linetype = "dashed") +
  geom_point( color=cols, size=2, alpha=0.9) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.title = element_text(size = 15,
                              color = "black",
                              face = "bold"),
    strip.text = element_text(face = "bold", color = "white", size = 12),
    strip.background = element_rect(fill = "black", linetype = "solid",
                                    color = "black", linewidth = 1)) +  
  facet_grid(num ~ ., scales = "free", space = "free") +
  labs(y = "Std. error")
  

#combined plot
plot_grid(asd, asd1, labels = "AUTO", ncol = 2,  rel_widths = c(4, 1.5), align="hv", hjust = 0, vjust = 1.5, scale = 1) #change names to a1, a2, a3 for each trial

#export summary tables
require(gt)

gt_tbl = res[,-1] #remove the abbreviated names
gt_tbl = gt_tbl[,c(7,6,5,3,4, 1, 2)] #reorder columns
colnames(gt_tbl)[colnames(gt_tbl) == 'num'] <- 'Missingness (percent)' #fix names
colnames(gt_tbl)[colnames(gt_tbl) == 'estimate'] <- 'Estimate'
colnames(gt_tbl)[colnames(gt_tbl) == 'std.error'] <- 'Std. error'
colnames(gt_tbl)[colnames(gt_tbl) == 'mech'] <- 'Missingness mechanism'
rownames(gt_tbl) = NULL #remove rownames
gt_tbl <- gt(gt_tbl)

gt_tab = gt_tbl |>  tab_header(
  title = md("**Simulation results**"),
  subtitle = md("Bias of the beta estimates, standard error and confidence intervals of different imputation results compared to the original full dataset")) |>
  cols_align(
    align = "center")
#output file name
gtsave(gt_tab, "table.docx") #change name to table 1, 2, 3 for each run

#after generating N plots combine them in a pdf of png or whatever

png("a4_output.png", paper="a4")
plot_grid(a1, a2, a3, labels = NULL, nrow = 3, align="hv")
dev.off()