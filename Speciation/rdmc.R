library(rdmc)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 15))
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/convergentEvo")
#load example data
data(neutral_freqs)
data(selected_freqs)
data(positions)

##这里是设置参数文件
param_list <-
  parameter_barge(
    Ne =  10000,
    rec = 0.005,
    neutral_freqs = neutral_freqs,
    selected_freqs = selected_freqs,
    selected_pops = c(1, 3, 5),
    positions = positions,
    n_sites = 10,
    sample_sizes = rep(10, 6),
    num_bins = 1000,
    sels = c(
      1e-4,
      1e-3,
      0.01,
      seq(0.02, 0.14, by = 0.01),
      seq(0.15, 0.3, by = 0.05),
      seq(0.4, 0.6, by = 0.1)
    ),
    times = c(0, 5, 25, 50, 100, 500, 1000, 1e4, 1e6),
    gs = c(1 / (2 * 10000), 10 ^ -(4:1)),
    migs = c(10 ^ -(seq(5, 1, by = -2)), 0.5, 1),
    sources = selected_pops,
    locus_name = "test_locus",
    cholesky = TRUE
  )

#fit composite likelihood models
neut_cle <- mode_cle(param_list, mode = "neutral")
ind_cle <- mode_cle(param_list, mode = "independent")
mig_cle <- mode_cle(param_list, mode = "migration")
#sv_cle <- mode_cle(param_list, mode = "standing_source")


#update barge to fit a  mixed-mode model
param_list <-
  update_mode(
    barge = param_list,
    sets = list(c(1, 3), 5),
    modes = c("standing_source", "independent"))

#fit mixed-mode model
#multi_svind <- mode_cle(param_list, "multi")

#update to another mixed-mode
param_list <-
  update_mode(
    barge = param_list,
    sets = list(c(1, 3), 5),
    modes =  c("migration", "independent"))

#fit mixed-mode model
multi_migind <- mode_cle(param_list, "multi")

#merge data frame of all fit models
all_mods <-
  bind_rows(
    ind_cle,
    mig_cle,
    #sv_cle,
    #multi_svind,
    multi_migind
  )

#max composite likelihood estimate 
#of all params over all models  
all_mods %>%
  group_by(model) %>%
  filter(cle == max(cle))

#write neutral and model results to separate files
readr::write_csv(neut_cle, "rdmc_neutral.csv")
readr::write_csv(all_mods, "rdmc_modes.csv")

neut <- unique(neut_cle$cle)
all_mods %>%
  group_by(selected_sites, model) %>%
  summarise(mcle = max(cle) - neut) %>%
  ggplot(aes(selected_sites, mcle, colour = model)) +
  geom_line() +
  geom_point() +
  xlab("Position") +
  ylab("Composite likelihood") +
  theme(legend.position = "n") +
  scale_color_brewer(palette = "Set1")

#visualize likelihood surface wrt selection coefficients
all_mods %>%
  group_by(sels, model) %>%
  summarise(mcle = max(cle) - neut) %>%
  ggplot(aes(sels, mcle, colour = model)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 0.05, lty = 2) +
  ylab("Composite likelihood") +
  xlab("Selection coefficient") +
  scale_color_brewer(palette = "Set1")



# Libraries
library(ggplot2)
library(hrbrthemes)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/permutation/XPCLR_recom_rate/norm20M_10k")
var11 = read.table("Chr_AB_15.xpclr.txt",header = F)
var1 = var11[var11$V1 == "chr1A",4]
var22 = read.table("Chr_AB_45.xpclr.txt",header = F)
var2 = var22[var22$V1 == "chr1A",4]

var1 = var1[sample(1:length(var1),1000)]
var2 = var2[sample(1:length(var2),1000)]
#hist(var1,breaks = 50)
# Dummy data
data <- data.frame(var1 ,var2)
# Chart
p <- ggplot(data, aes(x=x) ) +
  # Top
  geom_density( aes(x = var1, y = ..density..), fill="#69b3a2" ) +
  geom_label( aes(x=4.5, y=0.25, label="variable1"), color="#69b3a2") +
  # Bottom
  geom_density( aes(x = var2, y = -..density..), fill= "#404080") +
  geom_label( aes(x=4.5, y=-0.25, label="variable2"), color="#404080") +
  theme_ipsum() +
  xlab("value of x")
p

cdfPlot(var1)
plot(ecdf(EuStockMarkets[,"DAX"]))
plot(ecdf(iris$Sepal.Length))
######
varA = var1[var11$V4>1]
varB = var2[var11$V4>1]
varA = varA[sample(1:length(varA),1000)]
varB = varB[sample(1:length(varB),1000)]
#hist(var1,breaks = 50)
# Dummy data
data <- data.frame(varA ,varB)
# Chart
p <- ggplot(data, aes(x=x) ) +
  # Top
  geom_density( aes(x = varA, y = ..density..), fill="#69b3a2" ) +
  geom_label( aes(x=4.5, y=0.25, label="variable1"), color="#69b3a2") +
  # Bottom
  geom_density( aes(x = varB, y = -..density..), fill= "#404080") +
  geom_label( aes(x=4.5, y=-0.25, label="variable2"), color="#404080") +
  theme_ipsum() +
  xlab("value of x")
p











