#load libraries
library(tidyverse)

#import abundance values (major and subclass categories) and put in a dataframe
Major5LipidClasses<- read.csv('BP_input_mainv3.csv', header = TRUE, sep = ',')
Subclass22LipidClasses<- read.csv('BP_input_subclass.csv', header = TRUE, sep = ',')
NOFA_subclass<- read.csv('NOFA_subclass.csv', header = TRUE, sep = ',')
NOFA_Major<- read.csv('NOFA_Major.csv', header = TRUE, sep = ',')


Major5LipidClassesdf<-as.data.frame(Major5LipidClasses)
Subclass22LipidClassesdf<-as.data.frame(Subclass22LipidClasses)
NOFA_subclassdf<-as.data.frame(NOFA_subclass)
NOFA_Majordf<-as.data.frame(NOFA_Major)

#verify expected contents are present
summary(Major5LipidClassesdf)
summary(Subclass22LipidClassesdf)

#Closest to working. Need to fix axes
#Major Lipid Classes####
#By Raw abundance
ggplot(Major5LipidClassesdf, aes(x=Lipid_Main_Class, y=Raw.Abundance, color=Heat.Shock.Status))+ 
                            ylab("Raw Abundance")+ xlab("Lipid Class") +
                            geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.1)) +
                            theme(panel.background = element_rect(fill = "transparent")) +
                            ggtitle('Raw Abundance by Major Lipid Class') +
                            theme(plot.title = element_text(hjust = 0.5))

#BY PROPORTION
ggplot(Major5LipidClassesdf, aes(x=Lipid_Main_Class, y=Proportional.Abundance, color=Heat.Shock.Status))+ 
  ylab("Proportional Abundance")+ xlab("Lipid Class") +
  geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.1)) +
  theme(panel.background = element_rect(fill = "transparent")) +
  ggtitle('Proportional Abundance by Major Lipid Class') +
  theme(plot.title = element_text(hjust = 0.5))

#Subclass Lipids ####
#By Raw abundance
ggplot(Subclass22LipidClassesdf, aes(x=Lipid_Subclass, y=Raw.Abundance, color=Heat.Shock.Status))+ 
  ylab("Raw Abundance")+ xlab("Lipid Class") +
  geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.1)) +
  theme(panel.background = element_rect(fill = "transparent")) +
  ggtitle('Raw Abundance by Lipid Subclass') +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#By Raw abundance without NOFA
ggplot(NOFA_subclassdf, aes(x=Lipid_Subclass, y=Raw.Abundance, color=Heat.Shock.Status))+ 
  ylab("Raw Abundance")+ xlab("Lipid Class") +
  geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.1)) +
  theme(panel.background = element_rect(fill = "transparent")) +
  ggtitle('Raw Abundance by Lipid Subclass Without Fatty Acids') +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#BY PROPORTION
ggplot(Subclass22LipidClassesdf, aes(x=Lipid_Subclass, y=Proportional.Abundance, color=Heat.Shock.Status))+ 
  ylab("Proportional Abundance")+ xlab("Lipid Class") +
  geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.1)) +
  theme(panel.background = element_rect(fill = "transparent")) +
  ggtitle('Proportional Abundance by Lipid Subclass') +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#BY PROPORTION without PC
ggplot(NOFA_subclass, aes(x=Lipid_Subclass, y=Proportional.Abundance, color=Heat.Shock.Status))+ 
  ylab("Proportional Abundance")+ xlab("Lipid Class") +
  geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.1)) +
  theme(panel.background = element_rect(fill = "transparent")) +
  ggtitle('Proportional Abundance by Lipid Subclass Except Fatty Acids') +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Make changes here#####
ggplot(Major5LipidClassesdf, aes(x=Lipid_Main_Class, y=Raw.Abundance, color=Heat.Shock.Status))+ 
  ylab("Raw Abundance")+ xlab("Lipid Class") +
  geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.1)) +
  theme(panel.background = element_rect(fill = "transparent")) +
  ggtitle('Raw Abundance by Major Lipid Class') +
  theme(plot.title = element_text(hjust = 0.5))
