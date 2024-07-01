library(haven)
library(ggplot2)
library(ggalluvial)

#import the data
want <- read_sas("want.sas7bdat", NULL)

#drop the missing values
want_nomiss=subset(want,want$Speech2!='')

#select only the first four measurements
want_sub=subset(want_nomiss,want_nomiss$count<5)

#visualise the results
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


