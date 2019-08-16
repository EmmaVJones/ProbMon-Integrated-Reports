## Script to explore Cramer Von Mises test to compare 2 cdf curves

#install.packages('CDFt')
library(CDFt)

# Example script from https://www.r-bloggers.com/goodness-of-fit-test-in-r/
num_of_samples = 1000
x <- rgamma(num_of_samples, shape = 10, scale = 3)
x <- x + rnorm(length(x), mean=0, sd = .1)

num_of_samples = 100000
y <- rgamma(num_of_samples, shape = 10, scale = 3)
res <- CramerVonMisesTwoSamples(x,y)
p_value = 1/6*exp(-res)




# Test with our data 
library(plyr)
library(dplyr)
library(ggplot2)

# Bring in probmon_trend_status.csv
dat <- read.csv('data/probmon_trend_status.csv')


Bay_0107 <- filter(dat,Subpopulation=='Bay Watersheds 2001-2007'&Indicator=='VSCIAll')%>%select(Value,Estimate.P)
Bay_0814 <- filter(dat,Subpopulation=='Bay Watersheds 2008-2014'&Indicator=='VSCIAll')%>%select(Value,Estimate.P)

# Bay 01-07 vs post implementation
test <- CramerVonMisesTwoSamples(Bay_0107$Value,Bay_0814$Value)
1/6*exp(-test) #pvalue, not significant
ks.test(Bay_0107$Value,Bay_0814$Value) #pvalue significant


# CDF plot 
ggplot(Bay_0107,aes(Value,Estimate.P))+ geom_point()+
  labs(x='VSCI',y="Percentile") +
  ggtitle('Bay 2001-2007, 2008-2014') + 
  theme(plot.title = element_text(hjust=0.5,face='bold',size=15)) +
  theme(axis.title = element_text(face='bold',size=12))+
  geom_point(data=Bay_0107,aes(Value,Estimate.P),color='orange',size=2) + 
  geom_text(data=Bay_0107[50,],label='Bay 01-07',hjust=1.2,color='orange') +
  geom_point(data=Bay_0814,aes(Value,Estimate.P),color='gray',size=2)+
  geom_text(data=Bay_0814[75,],label='Bay 08-14',hjust=1.2,color='gray')



# nonbay different
nonBay_0107 <- filter(dat,Subpopulation=='Non-Bay Watersheds 2001-2007'&Indicator=='VSCIAll')%>%select(Value,Estimate.P)
nonBay_0814 <- filter(dat,Subpopulation=='Non-Bay Watersheds 2008-2014'&Indicator=='VSCIAll')%>%select(Value,Estimate.P)

test2 <- CramerVonMisesTwoSamples(nonBay_0107$Value,nonBay_0814$Value)
1/6*exp(-test2) #pvalue, not significant
ks.test(nonBay_0107$Value,nonBay_0814$Value)  #pvalue, not significant


# CDF plot 
ggplot(nonBay_0107,aes(Value,Estimate.P))+ geom_point()+
  labs(x='VSCI',y="Percentile") +
  ggtitle('nonBay 2001-2007, 2008-2014') + 
  theme(plot.title = element_text(hjust=0.5,face='bold',size=15)) +
  theme(axis.title = element_text(face='bold',size=12))+
  geom_point(data=nonBay_0107,aes(Value,Estimate.P),color='orange',size=2) + 
  geom_text(data=nonBay_0107[50,],label='nonBay 01-07',hjust=1.2,color='orange') +
  geom_point(data=nonBay_0814,aes(Value,Estimate.P),color='gray',size=2)+ 
  geom_text(data=nonBay_0814[75,],label='nonBay 08-14',hjust=1.2,color='gray')


# Bay nonbay 01-07 diff?
test3 <- CramerVonMisesTwoSamples(nonBay_0107$Value,Bay_0107$Value)
1/6*exp(-test3) #pvalue, not significant
ks.test(nonBay_0107$Value,Bay_0107$Value)#pvalue, not significant

# Bay nonbay 08-14 diff?
test4 <- CramerVonMisesTwoSamples(nonBay_0814$Value,Bay_0814$Value)
1/6*exp(-test4) #pvalue, not significant
ks.test(nonBay_0814$Value,Bay_0814$Value) #pvalue, not significant

# cramer is pretty stringent, test if cramer picks up any differences in cdf 
bigsandy <-  filter(dat,Subpopulation=='Big Sandy'&Indicator=='VSCIAll')%>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
cp <-  filter(dat,Subpopulation=='Clinch-Powell'&Indicator=='VSCIAll')%>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
holston <-  filter(dat,Subpopulation=='Holston'&Indicator=='VSCIAll')%>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)

test5 <- CramerVonMisesTwoSamples(bigsandy$Value,holston$Value)
1/6*exp(-test5) #pvalue, FIRST SIGNIFICANT
ks.test(bigsandy$Value,holston$Value)#pvalue, significant

# CDF plot 
ggplot(bigsandy,aes(Value,Estimate.P))+ geom_point()+
  labs(x='VSCI',y="Percentile") +
  ggtitle('big sandy vs holston VSCI') + 
  theme(plot.title = element_text(hjust=0.5,face='bold',size=15)) +
  theme(axis.title = element_text(face='bold',size=12))+
  geom_point(data=bigsandy,aes(Value,Estimate.P),color='orange',size=2) + 
  geom_text(data=bigsandy[10,],label='big sandy',hjust=1.2,color='orange') +
  geom_point(data=holston,aes(Value,Estimate.P),color='gray',size=2)+ 
  geom_text(data=holston[25,],label='holston',hjust=1.2,color='gray')


test6 <- CramerVonMisesTwoSamples(bigsandy$Value,cp$Value)
1/6*exp(-test6) #pvalue, not significant
ks.test(bigsandy$Value,cp$Value)#pvalue, not significant

# test does order matter?
test7 <-  CramerVonMisesTwoSamples(cp$Value,bigsandy$Value)
1/6*exp(-test7) #nope, same answer
ks.test(cp$Value,bigsandy$Value) #nope, same answer


