library("ggplot2")
library("magrittr")

# replicate Leimar & McNamara's Figure 4a

# read in the data
the.data <- read.table("../../data/summary_vary_qmat_vs_qjuv.csv"
                       ,sep=";"
                       ,header=T)

str(the.data)

# calculate total variance and maternal variance
the.data[,"maternal_variance"] <- with(the.data, var_component_amat_envt + var_component_amat_phen + cov_amat_envt_ajuv + cov_amat_phen_ajuv)
the.data[,"maternal_envtal_variance"] <- with(the.data, var_component_amat_envt + cov_amat_envt_ajuv)
the.data[,"total_variance"] <- the.data[,"maternal_variance"] + with(the.data, var_ajuv + cov_amat_envt_ajuv + cov_amat_phen_ajuv)

the.data[,"maternal_variance"] <- with(the.data, maternal_variance/total_variance)
the.data[,"maternal_envtal_variance"] <- with(the.data, maternal_envtal_variance/total_variance)
the.data[,"juvenile_variance"] <- with(the.data, 1.0 - maternal_variance)

the.data.m01 <- the.data %>% subset(m == 0.1) 

# calculate maternal variance component at the logit scale
# which, accordingly to Leimar is given by the maternal variance
# plus the covariance between maternal and ajuv
ggplot(data=the.data.m01
       ,aes(x=qmat)
) + 
geom_point(aes(y = var_component_amat, colour="Total maternal variance")) +
geom_point(aes(y = var_component_amat_envt, colour="Maternal environmental variance")) +
geom_point(aes(y = var_component_amat_phen, colour="Maternal phenotypic variance")) +
geom_point(aes(y = var_component_ajuv, colour="Juvenile cue variance")) +
geom_point(aes(y = var_component_gen, colour="Genetic cue variance")) +
xlab("Accuracy of maternal cue") +
ylab("Absolute variance") +
theme_classic()
ggsave("maternal_variance_vs_juvenile_variance.pdf")


ggplot(data=the.data.m01
       ,aes(x=qmat)
) + 
geom_point(aes(y = cov_amat_envt_ajuv, colour="Cov(Juvenile cue, mat envt)")) +
geom_point(aes(y = cov_amat_phen_ajuv, colour="Cov(Juvenile cue, mat phen)")) +
geom_point(aes(y = cov_amat_ajuv, colour="Cov(Juvenile cue, total mat)")) +
xlab("Accuracy of maternal cue") +
ylab("Absolute covariance") +
theme_classic()
ggsave("maternal_covariance_vs_juvenile_variance.pdf")

