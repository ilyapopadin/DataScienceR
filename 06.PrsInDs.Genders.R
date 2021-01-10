#### not yet FINAL !!!

rm(list=ls(all=TRUE))

AllTraits = read.table('../../body/2derived/TraitsTableGeneATLAS.csv', head = TRUE, sep='\t')

load('../../body/2derived/04.HarvesterPrsInDs.AllPrs.RData')
prs=Final

VecOfTraits = unique(prs$Phenotype); length(VecOfTraits)

OneLine=data.frame()
for (i in 1 : length(VecOfTraits)) 
{ # i = 1
  trait = prs[prs$Phenotype == VecOfTraits[i],]
  females = trait[trait$OffsprSex == 2,]; nrow(females)
  males = trait[trait$OffsprSex == 1,]; 
  males = males[sample(nrow(males), 288),]
  
  DaughtersDistanceToMothersVsFathers.4e.08 = median(abs(females$Prs.ObsMinMother.Pt_4e.08)) - median(abs(females$Prs.ObsMinFather.Pt_4e.08))
  DaughtersDistanceToMothersVsFathers.0.0005 = median(abs(females$Prs.ObsMinMother.Pt_0.0005)) - median(abs(females$Prs.ObsMinFather.Pt_0.0005))
  
  SonsDistanceToMothersVsFathers.4e.08 = median(abs(males$Prs.ObsMinMother.Pt_4e.08)) - median(abs(males$Prs.ObsMinFather.Pt_4e.08))
  SonsDistanceToMothersVsFathers.0.0005 = median(abs(males$Prs.ObsMinMother.Pt_0.0005)) - median(abs(males$Prs.ObsMinFather.Pt_0.0005))
  
  OneLine = rbind(OneLine,data.frame(VecOfTraits[i],DaughtersDistanceToMothersVsFathers.4e.08,DaughtersDistanceToMothersVsFathers.0.0005,SonsDistanceToMothersVsFathers.4e.08,SonsDistanceToMothersVsFathers.0.0005))
}

summary(OneLine$DaughtersDistanceToMothersVsFathers.0.0005)
summary(OneLine$DaughtersDistanceToMothersVsFathers.4e.08)
summary(OneLine$SonsDistanceToMothersVsFathers.0.0005)
summary(OneLine$SonsDistanceToMothersVsFathers.4e.08)

sd(OneLine$DaughtersDistanceToMothersVsFathers.0.0005)
sd(OneLine$SonsDistanceToMothersVsFathers.0.0005)

sd(OneLine$DaughtersDistanceToMothersVsFathers.4e.08)
sd(OneLine$SonsDistanceToMothersVsFathers.4e.08)

OneLine = merge(AllTraits,OneLine)


rm(list=ls(all=TRUE))
Trios = read.table('../../body/2derived/FULL_DS_TRIOS.txt', head = TRUE, sep = '\t')
table(Trios$OffsprSex) # 340 males & 288 females
Traits = read.table('../../body/2derived/TraitsTableGeneATLAS.csv', head = TRUE, sep='\t')

OneLine=c()
Dirs <- dir("../../body/2derived/PRSiceOutputs/")
for (i in 1 : length(Dirs)) 
{ # i = 1
  Path=paste('../../body/2derived/PRSiceOutputs/',Dirs[i],'/',Dirs[i],'.all_score',sep='')
  if (file.exists(Path))
  {
  Prs = read.table(Path, head = TRUE)
  Prs = Prs[c(2,3,4)]
  Prs1 = Prs
  # calculate deviation
  names(Prs1)=c(paste(names(Prs),'Offspr',sep='.'))
  Trio1 = merge(Trios, Prs1, by.x = 'OffsprID', by.y = 'IID.Offspr', all.x = TRUE)
  names(Prs1)=c(paste(names(Prs),'Mother',sep='.'))
  Trio2 = merge(Trio1, Prs1, by.x = 'MatID', by.y = 'IID.Mother', all.x = TRUE)
  names(Prs1)=c(paste(names(Prs),'Father',sep='.'))
  Trio3 = merge(Trio2, Prs1, by.x = 'PatID', by.y = 'IID.Father', all.x = TRUE)
  names(Trio3)
  # comparison of deviation from expectation for males and females without normalization: if selection truncates males more we expect the ratio RatioBetweenGendersNoNormal > 1 (females deviate more ~ have higher sd)
  Trio3$ObsMinExp.Pt_4e.08.NoNormalization = Trio3$Pt_4e.08.Offspr - (Trio3$Pt_4e.08.Mother + Trio3$Pt_4e.08.Father)/2;
  Trio3$ObsMinExp.Pt_0.0005.NoNormalization = Trio3$Pt_0.0005.Offspr - (Trio3$Pt_0.0005.Mother + Trio3$Pt_0.0005.Father)/2;
  table(Trio3$OffsprSex) # 340 males and 288 females
  
  #A = wilcox.test(sample(abs(Trio3[Trio3$OffsprSex == 1,]$ObsMinExp.Pt_0.0005.NoNormalization),288),abs(Trio3[Trio3$OffsprSex == 2,]$ObsMinExp.Pt_0.0005.NoNormalization)); PDiffBetweenGendersNoNormal = A$p.value
  RatioBetweenGendersNoNormal = median(abs(Trio3[Trio3$OffsprSex == 2,]$ObsMinExp.Pt_0.0005.NoNormalization))/median(sample(abs(Trio3[Trio3$OffsprSex == 1,]$ObsMinExp.Pt_0.0005.NoNormalization),288))
  RatioSDBetweenGendersNoNormal = median(sd(Trio3[Trio3$OffsprSex == 2,]$ObsMinExp.Pt_0.0005.NoNormalization)/sd(sample(Trio3[Trio3$OffsprSex == 1,]$ObsMinExp.Pt_0.0005.NoNormalization),288))
  
  Trio3$ObsMinExp.Pt_4e.08 = (Trio3$Pt_4e.08.Offspr - (Trio3$Pt_4e.08.Mother + Trio3$Pt_4e.08.Father)/2) / sd((Trio3$Pt_4e.08.Mother + Trio3$Pt_4e.08.Father)/2); summary(Trio3$ObsMinExp.Pt_4e.08)
  Trio3$ObsMinExp.Pt_0.0005 = (Trio3$Pt_0.0005.Offspr - (Trio3$Pt_0.0005.Mother + Trio3$Pt_0.0005.Father)/2) / sd((Trio3$Pt_0.0005.Mother + Trio3$Pt_0.0005.Father)/2); summary(Trio3$ObsMinExp.Pt_0.0005)
  Trio3=Trio3[c(3,6,15,16)]
  # Trio3 = Trio3[Trio3$OffsprSex == 1,] #### GENDER!!!
  # Trio3 = Trio3[!Trio3$OffsprID %in% OftenProblematicOffspring,]
  Trio3$Phenotype = Dirs[i]
  A= wilcox.test(Trio3[Trio3$OffsprSex == 1,]$ObsMinExp.Pt_4e.08,Trio3[Trio3$OffsprSex == 2,]$ObsMinExp.Pt_4e.08, alternative = 'greater'); PDiffBetweenGenders = A$p.value # '1' = male, '2' = female,
  MeanAbsDiffMales = mean(sample(abs(Trio3[Trio3$OffsprSex == 1,]$ObsMinExp.Pt_4e.08),288)) # the same number as girls
  MeanAbsDiffFemales = mean(abs(Trio3[Trio3$OffsprSex == 2,]$ObsMinExp.Pt_4e.08))
  A= wilcox.test(abs(Trio3[Trio3$OffsprSex == 1,]$ObsMinExp.Pt_4e.08),abs(Trio3[Trio3$OffsprSex == 2,]$ObsMinExp.Pt_4e.08)); PAbsDiffBetweenGenders = A$p.value
  A = t.test(Trio3$ObsMinExp.Pt_4e.08, mu = 0); TTestPvaluePt_4e.08 = A$p.value
  A = wilcox.test(Trio3$ObsMinExp.Pt_4e.08, mu = 0); WilcoxTestPvaluePt_4e.08 = A$p.value
  A = t.test(Trio3$ObsMinExp.Pt_0.0005, mu = 0); TTestPvaluePt_0.0005 = A$p.value
  A = wilcox.test(Trio3$ObsMinExp.Pt_0.0005, mu = 0); WilcoxTestPvaluePt_0.0005 = A$p.value
  OneLine = rbind(OneLine,c(Dirs[i],TTestPvaluePt_4e.08,TTestPvaluePt_0.0005,WilcoxTestPvaluePt_4e.08,WilcoxTestPvaluePt_0.0005,PDiffBetweenGenders,PAbsDiffBetweenGenders,MeanAbsDiffMales,MeanAbsDiffFemales,RatioBetweenGendersNoNormal,PDiffBetweenGendersNoNormal,RatioSDBetweenGendersNoNormal))
  if (i == 1) {Final = Trio3}
  if (i  > 1) {Final = rbind(Final,Trio3)}
  }
}
OneLine= data.frame(OneLine);
names(OneLine)=c('Key','TTestPvaluePt_4e.08','TTestPvaluePt_0.0005','WilcoxTestPvaluePt_4e.08','WilcoxTestPvaluePt_0.0005','PDiffBetweenGenders','PAbsDiffBetweenGenders','MeanAbsDiffMales','MeanAbsDiffFemales','RatioBetweenGendersNoNormal','PDiffBetweenGendersNoNormal','RatioSDBetweenGendersNoNormal')
OneLine$RatioBetweenGendersNoNormal = as.numeric(OneLine$RatioBetweenGendersNoNormal)
OneLine$RatioSDBetweenGendersNoNormal= as.numeric(OneLine$RatioSDBetweenGendersNoNormal)
summary(OneLine$RatioBetweenGendersNoNormal) # a bit more than one
summary(OneLine[OneLine$PDiffBetweenGendersNoNormal <0.1,]$RatioBetweenGendersNoNormal) # more shifted
plot(OneLine$RatioBetweenGendersNoNormal,OneLine$PDiffBetweenGendersNoNormal)
summary(log2(OneLine$RatioSDBetweenGendersNoNormal)) # a bit more than one, median = 0.2355
hist(log2(OneLine$RatioSDBetweenGendersNoNormal), breaks=100, main = "log2(sd(Obs-Exp Females)/sd(Obs - Exp Males)), N = 691 traits"); abline(v=0,col = 'blue', lwd = 3); abline(v=median(log2(OneLine$RatioSDBetweenGendersNoNormal)),col = 'red', lwd = 3);
wilcox.test(OneLine$RatioSDBetweenGendersNoNormal, mu = 0)



Agg = aggregate(list(Final$ObsMinExp.Pt_4e.08,Final$ObsMinExp.Pt_0.0005), by = list(Final$Phenotype), FUN = mean)
names(Agg)=c('Key','Median.ObsMinExp.Pt_4e.08','Median.ObsMinExp.Pt_0.0005')
Agg = merge(Agg, Traits, all.x=TRUE)
Agg = merge(Agg,OneLine, all.x = TRUE)
Agg=Agg[order(abs(Agg$Median.ObsMinExp.Pt_4e.08), decreasing = TRUE),]
Agg=Agg[order(Agg$WilcoxTestPvaluePt_0.0005),] # head(Agg)
Agg=Agg[order(Agg$TTestPvaluePt_0.0005),]
Agg=Agg[order(Agg$WilcoxTestPvaluePt_4e.08),]
Agg=Agg[order(Agg$PAbsDiffBetweenGenders),]
Agg=Agg[order(abs(Agg$Median.ObsMinExp.Pt_0.0005), decreasing = TRUE),]
names(Agg)
Agg1=Agg[,c(1,2,3,4,6,10,11)]

str(Agg)
Agg$FemalesToMales = as.numeric(Agg$MeanAbsDiffFemales)/as.numeric(Agg$MeanAbsDiffMales); 
summary(Agg$FemalesToMales)
hist(Agg$FemalesToMales, breaks = 100)
hist(log2(Agg$FemalesToMales), breaks = 100)
wilcox.test(Agg$FemalesToMales, mu = 1)
plot(Agg$FemalesToMales,Agg$PAbsDiffBetweenGenders)
Agg=Agg[order(Agg$FemalesToMales, decreasing = TRUE),]

## the minimal p-values (should be 0.05/700 = 7.142857e-05) 0.05/100 
min(Agg$WilcoxTestPvaluePt_0.0005) # 0.001
min(Agg$WilcoxTestPvaluePt_4e.08) # 0.0001
min(Agg$TTestPvaluePt_4e.08) # 0.00424
min(Agg$TTestPvaluePt_0.0005) # 0.00294


pdf('temp.pdf')
par(mfrow=c(2,1))
hist(Agg$Median.ObsMinExp.Pt_4e.08, breaks = 100)
hist(Agg$Median.ObsMinExp.Pt_0.0005, breaks = 3000, xlim = c(-1e-04,1e-04), main = 'median(deviation) for each trait')

# 23105-0.0 basal metabolic rate
# 20023-0.0 Mean time to correctly identify matches
# 23106-0.0 Impedance of whole body
# 30200-0.0 Neutrophill percentage
par(mfrow=c(2,1))
hist(Final[Final$Phenotype == '20023-0.0',]$ObsMinExp.Pt_4e.08, breaks = 100)
hist(Final[Final$Phenotype == '20023-0.0',]$ObsMinExp.Pt_0.0005, breaks = 100)
hist(Final[Final$Phenotype == '924-0.0',]$ObsMinExp.Pt_4e.08, breaks = 100)
hist(Final[Final$Phenotype == '1408-0.0',]$ObsMinExp.Pt_0.0005, breaks = 200)

dev.off()

##### IF CONTRASTS ARE SIMILAR WITH DIFFERENT THRESHOLDS? YES!!!
cor.test(Final$ObsMinExp.Pt_4e.08,Final$ObsMinExp.Pt_0.0005, method = 'spearman')
plot(Final$ObsMinExp.Pt_4e.08,Final$ObsMinExp.Pt_0.0005)
cor.test(Agg$Median.ObsMinExp.Pt_4e.08,Agg$Median.ObsMinExp.Pt_0.0005, method = 'spearman')
plot(Agg$Median.ObsMinExp.Pt_4e.08,Agg$Median.ObsMinExp.Pt_0.0005)


##### IF OUTLAYERS ARE THE SAME IN DIFFERENT TRAITS - THERE ARE SOME SUSPICIOUS!
Extremes=c()
VectorOfKeys = unique(Final$Phenotype); length(VectorOfKeys)
for (i in 1:length(VectorOfKeys))
{ # i = 1
  temp = Final[Final$Phenotype == VectorOfKeys[i],]
  temp = temp[temp$ObsMinExp.Pt_4e.08 < quantile(temp$ObsMinExp.Pt_4e.08, 0.01) | temp$ObsMinExp.Pt_4e.08 > quantile(temp$ObsMinExp.Pt_4e.08, 0.99),]
  Extremes=c(Extremes,temp$OffsprID)
}
Out = data.frame(table(Extremes)); Out=Out[order(Out$Freq, decreasing = TRUE),]
hist(Out$Freq, breaks = 100)

OftenProblematicOffspring = Out[Out$Freq > 100,]$Extremes
length(OftenProblematicOffspring) # 18

### BOYS AND GIRLS SEPARATELY!?

