library(tidyverse)
library(lubridate)
library(reshape2)
library(brms)
library(RColorBrewer)
library(patchwork)
library(ggbeeswarm)
library(viridis)
library(survival)
library(survminer)
library(rstanarm)
library(coxme)
library(lmtest)
theme_set(theme_classic(base_size = 16))

## LOAD DATA
setwd("~/Desktop/SPRF/Geladas/Infant Survival/GitHub")
offspring.survival<-read.csv("Feder et al_Offspring survival.csv")
grooming.by.focal<-read.csv("Feder et al_By focal.csv")
male.events<-read.csv("Feder et al_FM_grooming.csv")
fem.events<-read.csv("Feder et al_FF_grooming.csv")

offspring.survival$mom
model1<-lm(Index.F~scale(Total.Mom.Focal), data=offspring.survival)
model2<-lm(Index.M~scale(Total.Mom.Focal), data=offspring.survival)
bptest(model1)
bptest(model2)

## PLOT S1
## ADJUST SIZE OF POINTS TO MATERNAL FOCAL EFFORT
sizes<-sqrt(offspring.survival$Total.Mom.Focal)
plot1a<-ggplot(data=offspring.survival, aes(x=N.Females, y=Index.F.original)) + 
  geom_point(color="#1F78B4", alpha=0.3, size=sizes)  + geom_smooth(method="lm", color="black") +
  ylab("FF index - original") + xlab("# of females") + scale_y_continuous(limits = c(-2.5,3.5), breaks=seq(-3,3,1)) +
  scale_x_continuous(limits = c(0,12), breaks=seq(0,12,3)) 
plot1a

plot1b<-ggplot(data=offspring.survival, aes(x=N.Females, y=Index.M.original)) + 
  geom_point(color="#33A02C", alpha=0.3, size=sizes) + geom_smooth(method="lm", color="black") +
  ylab("") + xlab("# of females") + scale_y_continuous(limits = c(-2.5,3.5), breaks=seq(-3,3,1)) +
  scale_x_continuous(limits = c(0,12), breaks=seq(0,12,3)) 
plot1b

plot1c<-ggplot(data=offspring.survival, aes(x=N.Females, y=Index.F)) + 
  geom_point(color="#1F78B4", alpha=0.3, size=sizes)  + geom_smooth(method="lm", color="black") +
  ylab("FF index - adjusted") + xlab("# of females") + scale_y_continuous(limits = c(-2.5,3.5), breaks=seq(-3,3,1)) +
  scale_x_continuous(limits = c(0,12), breaks=seq(0,12,3)) 
plot1c

plot1d<-ggplot(data=offspring.survival, aes(x=N.Females, y=Index.M)) + 
  geom_point(color="#33A02C", alpha=0.3, size=sizes) + geom_smooth(method="lm", color="black") +
  ylab("") + xlab("# of females") + scale_y_continuous(limits = c(-2.5,3.5), breaks=seq(-3,3,1)) +
  scale_x_continuous(limits = c(0,12), breaks=seq(0,12,3)) 
plot1d

plot1a+plot1b+
  plot1c+plot1d + plot_annotation(tag_levels = "a")

## RUN BAYESIAN SURVIVAL MODELS
surv_object<-Surv(time=offspring.survival$End.Age, event=offspring.survival$Censored)

## FIT NET EFFECT MODEL
fit.noTO <- stan_surv(surv_object ~   Takeover + scale(Mom.Age) + Index.M + Index.F + (1|Mom) + (1|Unit) + (1|Year), 
                      data = offspring.survival, cores=2, basehaz = 'weibull')
summary(fit.noTO,prob=c(0.055, 0.945), digits=3)

## FIT INTX EFFECT MODEL
fit <- stan_surv(surv_object ~   tve(Takeover, degree=3) + tve(Takeover, degree=3)*Index.F + tve(Takeover, degree=3)*Index.M + scale(Mom.Age) +  (1|Mom) + (1|Unit) + (1|Year),  
                 data = offspring.survival, cores=2, basehaz = 'weibull')
summary(fit, prob=c(0.055, 0.945), digits=2)
plot(fit, plotfun="tve",pars="TakeoverY")

## PLOT MODEL PREDICTIONS
## CREATE COLOR SCHEME
colors<-brewer.pal(4, "Paired")
male<-colors[c(4,3)]
females<-colors[c(2,1)]

## CREATE COMPLEMENTARY COLOR SCHEMES
female.colors<-RColorBrewer::brewer.pal(5,"Blues")
male.colors<-RColorBrewer::brewer.pal(5,"Greens")
female.colors<-female.colors[3:5]
male.colors<-male.colors[3:5]

## PLOT FEMALE-MALE EFFECTS
theme_set(theme_classic(base_size = 16))
quantile(offspring.survival$Index.M, probs=c(0.25, 0.75))

offspring.survival1 <- dplyr::select(offspring.survival, Mom.Age, Takeover, Index.M, Index.F, Unit, Year, Mom) %>%
  tidyr::expand(Index.M=c(-0.62,0.44), Takeover=c("N","Y"),
                nesting(Mom.Age, Index.F, Unit, Year, Mom)) %>%
  dplyr::mutate(id = 1:n(),
                Above.Male2 = as.factor(Index.M),
                Takeover2=as.factor(Takeover)) %>%
  nest_legacy(-c(Above.Male2, Takeover2))  %>%
  dplyr::mutate(., preds = map(
    data,
    ~ posterior_survfit(
      fit, prob=0.89,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))


offspring.survival1 <- unnest(offspring.survival1, preds) %>%
  dplyr::select(-data)

levels(offspring.survival1$Above.Male2)

library(viridis)
viridis(3)
offspring.survival1$Above.Male2<-factor(offspring.survival1$Above.Male2, levels=c("0.44", "-0.62"))
offspring.survival1.TO<-offspring.survival1[offspring.survival1$Takeover2=="Y",]
offspring.survival1.NoTo<-offspring.survival1[offspring.survival1$Takeover2=="N",]

# infants.social$Index.M2<-factor(infant.fissions$Index.M2, levels=c("","N"))

plot1c<-ggplot(offspring.survival1.NoTo,
               aes(
                 time,
                 fill = Above.Male2,
                 group = Above.Male2
               )) + labs(fill="female-male grooming") + 
  scale_fill_manual(labels = c("strong", "weak"),values = male) +
  scale_color_manual(values = male.colors[3:2]) +
  theme(legend.position = c(0.44,0.2)) +
  geom_ribbon(
    aes(time, ymin = ci_lb, ymax = ci_ub),
    col = NA,
    alpha = 0.7
  ) + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold")) +
  geom_line(aes(x=time, y = median, color=Above.Male2), lwd=1) +
  guides(color=F) +
  ylim(c(0, 1)) + xlim(0,5.0) + 
  xlab("offspring age (years)") + ylab("proportion of offspring surviving") 

plot1d<-ggplot(offspring.survival1.TO,
               aes(
                 time,
                 fill = Above.Male2,
                 group = Above.Male2
               )) + labs(fill="female-male grooming") + theme(legend.position = "none") +
  scale_fill_manual(values = male) +
  scale_color_manual(values = male.colors[3:2]) +
  geom_ribbon(
    aes(time, ymin = ci_lb, ymax = ci_ub),
    col = NA,
    alpha = 0.7
  ) + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold")) +
  geom_line(aes(x=time, y = median, color=Above.Male2), lwd=1) +
  guides(color=F) +
  ylim(c(0, 1)) + xlim(0,5.0) + 
  xlab("offspring age (years)")  + ylab("") 

plot1c+plot1d

## PLOT FEMALE EFFECTS
quantile(offspring.survival$Index.F, probs=c(0.25, 0.75))

offspring.survival1 <- dplyr::select(offspring.survival, N.Females, Mom.Age, Takeover, Index.M, Index.F, Unit, Year, Mom) %>%
  tidyr::expand(Index.F=c(.58,-.77), Takeover=c("N","Y"),
                nesting(N.Females, Mom.Age, Index.M, Unit, Year, Mom)) %>%
  dplyr::mutate(id = 1:n(),
                Above.Female2 = as.factor(Index.F),
                Takeover2=as.factor(Takeover)) %>%
  nest_legacy(-c(Above.Female2, Takeover2)) %>%
  dplyr::mutate(., preds = map(
    data,
    ~ posterior_survfit(
      fit, prob=0.89,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))

offspring.survival1 <- unnest(offspring.survival1, preds) %>%
  dplyr::select(-data)

library(viridis)
viridis(3)

offspring.survival1$Above.Female2<-factor(offspring.survival1$Above.Female2, levels=c("0.58", "-0.77"))
offspring.survival1.TO<-offspring.survival1[offspring.survival1$Takeover2=="Y",]
offspring.survival1.NoTO<-offspring.survival1[offspring.survival1$Takeover2=="N",]

plot1a<-ggplot(offspring.survival1.NoTO,
               aes(
                 time,
                 fill = Above.Female2,
                 group = Above.Female2
               )) + labs(fill="female-female grooming") + 
  scale_fill_manual(labels = c("strong", "weak"),values = females) +
  scale_color_manual(values = female.colors[c(3,2)]) +
  theme(legend.position = c(0.47,0.2)) +
  geom_ribbon(
    aes(time, ymin = ci_lb, ymax = ci_ub),
    col = NA,
    alpha = 0.7
  ) +  
  geom_line(aes(x=time, y = median, color=Above.Female2),lwd=1) +
  guides(color=F) +
  ylim(c(0, 1)) + xlim(0,5.0) +labs(title="non-takeover infants") +
  xlab("") + ylab("proportion of offspring surviving") + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title = element_text(face="bold"))

plot1b<-ggplot(offspring.survival1.TO,
               aes(
                 time,
                 fill = Above.Female2,
                 group = Above.Female2
               )) + labs(fill="female-female grooming") + theme(legend.position = "none") +
  scale_fill_manual(values = females) +
  scale_color_manual(values = female.colors[c(3,2)]) +
  geom_ribbon(
    aes(time, ymin = ci_lb, ymax = ci_ub),
    col = NA,
    alpha = 0.7
  ) +  labs(title="takeover infants") + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold")) +
  geom_line(aes(x=time, y = median, color=Above.Female2), lwd=1) +
  guides(color=F) +
  ylim(c(0, 1)) + xlim(0,5.0) +
  xlab("") + ylab("")

plot1a+plot1b + plot1c + plot1d +plot_annotation(tag_levels = "a")
ggsave(file="Gelada_bond_infant_survival_DYADIC.jpg", units="cm", width=22, height=20, dpi=300)

##
ggplot(data=offspring.survival, aes(x=Total.Mom.Focal,y=Index.F)) + geom_point()
ggplot(data=offspring.survival, aes(x=Total.Mom.Focal,y=Index.M)) + geom_point()

## PLOT RAW SURVIVAL CURVES
offspring.survival$fem.quartile <- dplyr::ntile(offspring.survival$Index.F, 3)  
offspring.survival$male.quartile <- dplyr::ntile(offspring.survival$Index.M, 3)  

no.takeover.offspring<-offspring.survival[offspring.survival$Takeover=='N',]
takeover.offspring<-offspring.survival[offspring.survival$Takeover=='Y',]
surv_object1<-Surv(time=no.takeover.offspring$End.Age, event=no.takeover.offspring$Censored)
surv_object2<-Surv(time=takeover.offspring$End.Age, event=takeover.offspring$Censored)

## GENERATE SURVPLOTS
## PLOT 1
survplot1<-survfit(surv_object1~fem.quartile, data=no.takeover.offspring)
survplot1<-ggsurvplot(survplot1, ylim=c(0,1), palette=female.colors, xlab="", ylab="Proportion surviving", legend="none", censor=F, title="no takeover")
survplot1<-survplot1$plot + theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

## PLOT 2
survplot2<-survfit(surv_object2~fem.quartile, data=takeover.offspring)
survplot2<-ggsurvplot(survplot2, ylim=c(0,1), palette=female.colors, xlab="", ylab="", legend="right",censor=F, legend.labs=c("Low","Mid","High"), legend.title="FF strength", title="takeover")
survplot2<-survplot2$plot + theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

## PLOT 3
survplot3<-survfit(surv_object1~male.quartile, data=no.takeover.offspring)
survplot3<-ggsurvplot(survplot3, ylim=c(0,1), palette=male.colors,  xlab="Age (years)",ylab="Proportion surviving", legend="none",censor=F)
survplot3<-survplot3$plot + theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

## PLOT 4
survplot4<-survfit(surv_object2~male.quartile, data=takeover.offspring)
survplot4<-ggsurvplot(survplot4, ylim=c(0,1), palette=male.colors, xlab="Age (years)", ylab="",legend="right",censor=F, legend.labs=c("Low","Mid","High"), legend.title="FM strength")
survplot4<-survplot4$plot + theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

library(patchwork)
(survplot1+survplot2)/(survplot3 + survplot4) + plot_annotation(tag_levels = "a")
ggsave("ProcB_FigureS3.jpg", units="cm", dpi=300, width=20, height=15) 

## PLOT 3
survplot.all1<-survfit(surv_object~fem.quartile, data=offspring.survival)
survplot.all1<-ggsurvplot(survplot.all1, ylim=c(0,1), palette=female.colors,  xlab="Age (years)",ylab="Proportion surviving", legend="none",censor=F)
survplot.all1<-survplot.all1$plot + theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

## PLOT 4
survplot.all2<-survfit(surv_object~male.quartile, data=offspring.survival)
survplot.all2<-ggsurvplot(survplot.all2, ylim=c(0,1), palette=male.colors,  xlab="Age (years)",ylab="", legend="none",censor=F)
survplot.all2<-survplot.all2$plot + theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
survplot.all1+survplot.all2 + plot_annotation(tag_levels = "a")

## MODEL CHANGES IN GROOMING EFFORT
## TOTAL GROOMING
grooming.by.focal$Infant<-as.character(grooming.by.focal$Infant)
grooming.by.focal$Period<-factor(grooming.by.focal$Period, levels=c("Pre", "Post"))

hist(grooming.by.focal$duration)
nrow(grooming.by.focal[(grooming.by.focal$duration<16 & grooming.by.focal$duration>14),])/
nrow(grooming.by.focal)

colnames(grooming.by.focal)
equation1<-brmsformula(Prop.F|weights(duration/15) ~ Infant*Period + (1|individual_id) + (1|TakeoverID), family="zero_inflated_beta")

hist(grooming.by.focal$duration)
model1<-brm(data=grooming.by.focal,Groom.F ~ Infant*Period + offset(log(duration)) + (1+Period|individual_id) + (1|TakeoverID), 
            family="hurdle_gamma", prior = c(
              prior(normal(0, 5), "Intercept"),
              prior(normal(0, 5), "b")),
            cores=2, control = list(adapt_delta=0.99))
summary(model1, prob=0.89)
conditional_effects(model1, "Infant:Period",prob = 0.89)

model2<-brm(data=grooming.by.focal,Groom.M ~ Infant*Period + offset(log(duration)) + (1+Period|individual_id) + (1|TakeoverID), 
            family="hurdle_gamma",prior = c(
              prior(normal(0, 5), "Intercept"),
              prior(normal(0, 5), "b")),
            cores=2, control = list(adapt_delta=0.99))
summary(model2, prob=0.89)
conditional_effects(model2, "Infant:Period",prob = 0.89)

sum(grooming.by.focal$Groom.M[grooming.by.focal$Period=="Pre" & grooming.by.focal$Infant=="Y"])/sum(grooming.by.focal$duration[grooming.by.focal$Period=="Pre" & grooming.by.focal$Infant=="Y"])
sum(grooming.by.focal$Groom.M[grooming.by.focal$Period=="Post" & grooming.by.focal$Infant=="Y"])/sum(grooming.by.focal$duration[grooming.by.focal$Period=="Post" & grooming.by.focal$Infant=="Y"])
sum(grooming.by.focal$Groom.M[grooming.by.focal$Period=="Pre" & grooming.by.focal$Infant=="N"])/sum(grooming.by.focal$duration[grooming.by.focal$Period=="Pre" & grooming.by.focal$Infant=="N"])
sum(grooming.by.focal$Groom.M[grooming.by.focal$Period=="Post" & grooming.by.focal$Infant=="N"])/sum(grooming.by.focal$duration[grooming.by.focal$Period=="Post" & grooming.by.focal$Infant=="N"])

## MAKE PLOT
to.plota<-conditional_effects(model1, prob = 0.89)
to.plotb<-conditional_effects(model2, prob = 0.89)

to.plota$`Infant:Period`
to.plota$`Infant:Period`$estimate__<-to.plota$`Infant:Period`$estimate__/to.plota$`Infant:Period`$duration
to.plota$`Infant:Period`$lower__<-to.plota$`Infant:Period`$lower__/to.plota$`Infant:Period`$duration
to.plota$`Infant:Period`$upper__<-to.plota$`Infant:Period`$upper__/to.plota$`Infant:Period`$duration

to.plotb$`Infant:Period`
to.plotb$`Infant:Period`$estimate__<-to.plotb$`Infant:Period`$estimate__/to.plotb$`Infant:Period`$duration
to.plotb$`Infant:Period`$lower__<-to.plotb$`Infant:Period`$lower__/to.plotb$`Infant:Period`$duration
to.plotb$`Infant:Period`$upper__<-to.plotb$`Infant:Period`$upper__/to.plotb$`Infant:Period`$duration

to.plotb$Period
to.plotb$`Period`$estimate__<-to.plotb$`Period`$estimate__/to.plotb$`Period`$duration
to.plotb$`Period`$lower__<-to.plotb$`Period`$lower__/to.plotb$`Period`$duration
to.plotb$`Period`$upper__<-to.plotb$`Period`$upper__/to.plotb$`Period`$duration

colors<-brewer.pal(4, "Paired")
male<-colors[c(3,4)]
females<-colors[c(1,2)]

to.plota<-plot(to.plota, plot = FALSE)[[3]] + scale_shape_discrete() + scale_color_manual(values = females)  +
  ylab("proportion of\nfocal time grooming") + labs(title="with females") + xlab("infant present?") + scale_y_continuous(limits=c(0,0.11), labels = seq(0,0.1,0.05), breaks=seq(0,0.1,0.05)) +
  theme(plot.title=element_text(face="bold", hjust=0.5), axis.title = element_text(face="bold"), legend.title=element_text(face="bold")) + guides(fill = "none")
to.plotb<-plot(to.plotb, plot = FALSE)[[3]] + scale_color_manual(values = male)  + xlab("infant present?") +
  ylab("") + scale_y_continuous(limits = c(0,0.11), labels = seq(0,0.1,0.05), breaks=seq(0,0.1,0.05))+ labs(title="with males") + 
  theme(plot.title=element_text(face="bold", hjust=0.5), axis.title = element_text(face="bold"), legend.title=element_blank())  + guides(fill = "none")

to.plota+to.plotb + plot_annotation(tag_levels="a") + plot_layout(guides = "collect")

to.plota$data

## TRY DIFFERENT WAY
library(scales)

grooming.by.focal$Prop.F<-round(grooming.by.focal$Groom.F/grooming.by.focal$duration,3)
grooming.by.focal$Zero.F<-"Yes"
grooming.by.focal$Zero.F[grooming.by.focal$Prop.F>0]<-"No"
grooming.by.focal$Prop.F[grooming.by.focal$Prop.F==0]<-0.001
grooming.by.focal$Prop.M<-round(grooming.by.focal$Groom.M/grooming.by.focal$duration,3)
grooming.by.focal$Zero.M<-"Yes"
grooming.by.focal$Zero.M[grooming.by.focal$Prop.M>0]<-"No"
grooming.by.focal$Prop.M[grooming.by.focal$Prop.M==0]<-0.001

jitter <- position_jitter(width = 0.1, height = 0.3)

?geom_po
plot.to.diff1<-ggplot() + geom_jitter(data=grooming.by.focal, aes(x=Period, y=Prop.F, color=Infant, alpha=Zero.F), width=0.3, height=0.3) +
  scale_alpha_discrete(range=c(0.5,0.05)) +
  facet_grid(cols = vars(Infant), labeller=as_labeller(c("N"="no infant","Y"="infant"))) + 
  geom_point(data=to.plota$data, aes(x=Period, y=estimate__), cex=2.5) + labs(title="with females") +
  scale_y_log10(labels = function(x) sprintf("%.2f",x)) +  # Use scale_y_log10 for log transformation
  geom_errorbar(data=to.plota$data, aes(x=Period, ymin=lower__, ymax=upper__), width=0, linewidth=1.5) +
  theme(plot.title=element_text(face="bold", hjust=0.5), axis.title = element_text(face="bold"), legend.title=element_blank())+ guides(color = "none", alpha="none") +
  scale_shape_discrete() + scale_color_manual(values = females) + ylab("proportion of\nfocal time grooming")

plot.to.diff2<-ggplot() + geom_jitter(data=grooming.by.focal, aes(x=Period, y=Prop.M, color=Infant,alpha=Zero.M), width=0.3, height=0.3) +
  scale_alpha_discrete(range=c(0.5,0.05)) +
  facet_grid(cols = vars(Infant), labeller=as_labeller(c("N"="no infant","Y"="infant"))) + 
  geom_point(data=to.plotb$data, aes(x=Period, y=estimate__), cex=2.5) + labs(title="with males") +
  scale_y_log10(labels = function(x) sprintf("%.2f",x)) +  # Use scale_y_log10 for log transformation
  geom_errorbar(data=to.plotb$data, aes(x=Period, ymin=lower__, ymax=upper__), width=0, linewidth=1.5) +
  theme(plot.title=element_text(face="bold", hjust=0.5), axis.title = element_text(face="bold"), legend.title=element_blank())  + guides(color = "none", alpha="none") +
  scale_shape_discrete() + scale_color_manual(values = male) + ylab("proportion of\nfocal time grooming")
  
plot.to.diff1/plot.to.diff2 + plot_annotation(tag_levels = "a")
ggsave(file="PLOT_FOR_REVIEWERS.jpg", units="cm", width=20, height=25, dpi=300)

## MODEL GROOMING EVENTS
## FEMALE-FEMALE EVENTS
fem.events$Period<-factor(fem.events$Period, levels=c("Pre","Post"))
male.events$Period<-factor(male.events$Period, levels=c("Pre","Post"))

model3<-brm(data=fem.events, Risk|weights(Duration)~ 
              Period*Infant + (1|CodeFemale) + (1|focal_unit_name/TakeoverID), 
              prior = c(
              prior(normal(0, 2), "Intercept"),
              prior(normal(0, 5), "b")),
              family="bernoulli", cores=4, control = list(adapt_delta=0.99))
summary(model3, prob=0.89)

conditional_effects(model3, "Period",prob=0.89)

## FEMALE-MALE EVENTS
model4<-brm(data=male.events, Risk|weights(Duration)~ 
              Period*Infant + N.Males + (1|CodeFemale) + (1|focal_unit_name/TakeoverID), 
            prior = c(
              prior(normal(0, 2), "Intercept"),
              prior(normal(0, 5), "b")),
            family="bernoulli", cores=4, control = list(adapt_delta=0.99))
summary(model4, prob = 0.89)

conditional_effects(model4)

##
male.colors<-c("palegreen1", "palegreen3","palegreen4" )
female.colors<-c("#165A85","#238dd0","#A6CEE3","darkgray")

## SUMMARIZE FF EVENTS
groom.summary.FF<-fem.events %>%
  group_by(Period, Infant, Partner.Cat) %>%
  dplyr::summarise(Grooming.Time=sum(Duration, na.rm = T)/3600)
groom.summary.FF$Period<-factor(groom.summary.FF$Period, levels=c("Pre", "Post"))
groom.summary.FF$Partner.Cat<-factor(groom.summary.FF$Partner.Cat, levels=c("Mother-daughter","Sisters","Other close kin","Non-kin"))

groom.summary.FF$Prop<-NA
groom.summary.FF$Prop[1:4]<-groom.summary.FF$Grooming.Time[1:4]/sum(groom.summary.FF$Grooming.Time[1:4])
groom.summary.FF$Prop[5:8]<-groom.summary.FF$Grooming.Time[5:8]/sum(groom.summary.FF$Grooming.Time[5:8])
groom.summary.FF$Prop[9:12]<-groom.summary.FF$Grooming.Time[9:12]/sum(groom.summary.FF$Grooming.Time[9:12])
groom.summary.FF$Prop[13:16]<-groom.summary.FF$Grooming.Time[13:16]/sum(groom.summary.FF$Grooming.Time[13:16])

## SUMMARIZE FM EVENTS
groom.summary.FM<-male.events %>%
  group_by(Period, Infant,Partner.Cat) %>%
  dplyr::summarise(Grooming.Time=sum(Duration, na.rm = T)/3600)
groom.summary.FM$Period<-factor(groom.summary.FM$Period, levels=c("Pre", "Post"))

groom.summary.FM$Prop<-NA
groom.summary.FM$Prop[1:3]<-groom.summary.FM$Grooming.Time[1:3]/sum(groom.summary.FM$Grooming.Time[1:3])
groom.summary.FM$Prop[4:6]<-groom.summary.FM$Grooming.Time[4:6]/sum(groom.summary.FM$Grooming.Time[4:6])
groom.summary.FM$Prop[7:9]<-groom.summary.FM$Grooming.Time[7:9]/sum(groom.summary.FM$Grooming.Time[7:9])
groom.summary.FM$Prop[10:12]<-groom.summary.FM$Grooming.Time[10:12]/sum(groom.summary.FM$Grooming.Time[10:12])

## PLOTS
plot3a<-ggplot(groom.summary.FF[groom.summary.FF$Infant=="N",], aes(x=Period, y=Prop, fill=Partner.Cat)) + geom_bar(stat="identity", alpha=0.7) + scale_fill_manual(values=female.colors) + 
  ylab("proportion of\nFF grooming effort") + xlab("") + ggtitle("infant absent") +
  labs(fill="Partner category") + theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold"), plot.title=element_text(face="bold", hjust = 0.5)) +
  theme(legend.position = "none")

plot3b<-ggplot(groom.summary.FF[groom.summary.FF$Infant=="Y",], aes(x=Period, y=Prop, fill=Partner.Cat)) + geom_bar(stat="identity", alpha=0.7) + 
  scale_fill_manual(values=female.colors) + ylab("") + xlab("")+
  labs(fill="Partner category") + ggtitle("infant present") + 
  theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold"), 
        plot.title=element_text(face="bold",hjust = 0.5))

plot3c<-ggplot(groom.summary.FM[groom.summary.FM$Infant=="N",], aes(x=Period, y=Prop, fill=Partner.Cat)) + geom_bar(stat="identity", alpha=0.9) + scale_fill_manual(values=male.colors) + 
  ylab("proportion of\nFM grooming effort") +
  labs(fill="Partner category") + theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold")) +
  theme(legend.position = "none")

plot3d<-ggplot(groom.summary.FM[groom.summary.FM$Infant=="Y",], aes(x=Period, y=Prop, fill=Partner.Cat)) + geom_bar(stat="identity", alpha=0.9) + scale_fill_manual(values=male.colors) + ylab("") + 
  labs(fill="Partner category") + theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold"))

plot3a+plot3b +
  plot3c+plot3d + plot_annotation(tag_levels = "a")

## FIGURE S2
## VISUALIZED AT FOCAL LEVEL
## TRY OBSERVATION-LEVEL
fem.events$focal_unit_name
fem.events$Risk
groom.by.id<-fem.events %>%
  group_by(CodeFemale, Infant, Period, TakeoverID) %>%
  dplyr::summarise(Risk.Prop=sum(Duration[Risk==1], na.rm = T)/sum(Duration, na.rm = T), Duration=sum(Duration, na.rm = T))

groom.by.id$Point.Size<-sqrt(groom.by.id$Duration)

effect1<-as.data.frame(conditional_effects(model3, "Period:Infant", prob=0.89)$Period)

library(ggbeeswarm)
plots2a<-ggplot() + theme(legend.position = "none") +
  geom_quasirandom(data=groom.by.id[groom.by.id$Infant=="N",], aes(x=Period, y=Risk.Prop, size = Point.Size),
                   alpha=0.4, color=females[1])  + ylab("proportion of\nFF grooming w/ non-kin") +
  geom_point(data=effect1[effect1$Infant=="N",], aes(x=Period, y=estimate__))+
  geom_errorbar(data=effect1[effect1$Infant=="N",], aes(x=Period, ymin=lower__, ymax=upper__), width=0.1, lwd=1) +
  ggtitle("infant absent") +  theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold"), plot.title=element_text(face="bold",hjust = 0.5))+
  xlab("")

plots2b<-ggplot() + theme(legend.position = "none") +
  geom_quasirandom(data=groom.by.id[groom.by.id$Infant=="Y",], aes(x=Period, y=Risk.Prop, size = Point.Size),
                   alpha=0.4, color=females[2])  + ylab("") +
  geom_point(data=effect1[effect1$Infant=="Y",], aes(x=Period, y=estimate__))+
  geom_errorbar(data=effect1[effect1$Infant=="Y",], aes(x=Period, ymin=lower__, ymax=upper__), width=0.1, lwd=1) +
  ggtitle("infant present") +  theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold"), plot.title=element_text(face="bold",hjust = 0.5))+
  xlab("")

groom.by.id2<-male.events %>%
  group_by(CodeFemale, Infant, Period, TakeoverID) %>%
  dplyr::summarise(Risk.Prop=sum(Duration[Risk==1], na.rm = T)/sum(Duration, na.rm = T), Duration=sum(Duration, na.rm = T))

groom.by.id2$Point.Size<-sqrt(groom.by.id2$Duration)

effect2<-as.data.frame(conditional_effects(model4, "Period:Infant", prob=0.89)$Period)

plots2c<-ggplot() + theme(legend.position = "none") +
  geom_quasirandom(data=groom.by.id2[groom.by.id2$Infant=="N",], aes(x=Period, y=Risk.Prop, size = Point.Size),
                   alpha=0.4, color=male[1]) + 
  geom_point(data=effect2[effect2$Infant=="N",], aes(x=Period, y=estimate__))+
  geom_errorbar(data=effect2[effect2$Infant=="N",], aes(x=Period, ymin=lower__, ymax=upper__), width=0.1, lwd=1) +
  ylab("proportion of\nFM grooming w/ leader male") + theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold"), plot.title=element_text(face="bold",hjust = 0.5))

plots2d<-ggplot() + theme(legend.position = "none") +
  geom_quasirandom(data=groom.by.id2[groom.by.id2$Infant=="Y",], aes(x=Period, y=Risk.Prop, size = Point.Size),
                   alpha=0.4, color=male[2]) + 
  geom_point(data=effect2[effect2$Infant=="Y",], aes(x=Period, y=estimate__))+
  geom_errorbar(data=effect2[effect2$Infant=="Y",], aes(x=Period, ymin=lower__, ymax=upper__), width=0.1, lwd=1) +
  ylab("") + theme(axis.title = element_text(face="bold"), legend.title=element_text(face="bold"), plot.title=element_text(face="bold",hjust = 0.5))

(plots2a+plots2b)/(plots2c+plots2d) + plot_annotation(tag_levels = "a")
ggsave(file="Gelada_partner_choice.jpg", units="cm", width=22, height=20, dpi=300)

