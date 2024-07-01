library(haven)
library(ggplot2)
library(ggalluvial)
setwd(dir = 'C:/Users/u0118563/OneDrive - KU Leuven/Projecten/Rare diseases/Data')
want <- read_sas("want.sas7bdat", NULL)
want_nomiss=subset(want,want$Speech2!='')
want_sub=subset(want_nomiss,want_nomiss$count<5)
want_sub$Speech=want_sub$Speech2
ggplot(want_sub,
       aes(x = count, stratum = Speech, alluvium = id,
           fill = Speech, label = Speech)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  xlab('# Measurement')+ 
  theme(text = element_text(size=20),legend.title = element_blank(),legend.text=element_text(size=rel(1.1)))+
  theme_minimal()


