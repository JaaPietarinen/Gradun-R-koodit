#Liite 
#Ohjelmat genominlaajuisen kartoituksen tekemiseksi muunnetulla BGLR-paketilla. 
#Samat ohjelmat tulevat jakoon dokumentointeineen mahdollisessa artikkelissa, 
#joten dokumentaatio ja funktioiden kuvaukset ovat englanniksi.
#
#
#For easy readability consider opening this file (code) in Notepad++ (https://notepad-plus-plus.org/) and then changing the language from the settings into R
#
#In this documantation are provided the R-codes required to run GWAS with Sure independence screening (SIS)
#(Fan & Lv 2008) using the BGLR-package (Pèrez & de los Campos 2014)that has been modified to also calculate
# statistics for modified Bayes factors. Codes for producing random test
#datasests and simulated phenotypes are also provided as steps 0.1 and 0.2. In this pipeline the program SHAPEIT
#(Delaneu et. al. 2012) is also used for imputing missing data as step 3. in the pipeline. For SHAPEIT only an
#example of a command line is provided. Go to www.shapeit.fr/ for more information and manuals. 
#
#In this pipeline the format of genotype files has six columns of ID and pedigree information after which follow the SNP information columns, 
#where each cell contains both alleles for the SNP separated by space. All the files are tab delimited. 
#The map file has four columns where the first has information about chromosome position,
#the second contains the names of the SNPs, the third contains information about CM distance but can be zero and the fourth
#is a column of the nucleotide positions.
#
#A function to create purely random data set for testing the pipeline is provided in the beginning
#of the documentation (Step 0.1 and 0.2).
#
#Please report bugs or other issues to jaakko.pietarinen@gmail.com
#--------------------------------------------------------------------------------------------------------------

#Step 0.1
#
#Creating a purely random dataset (genotype and map) for testing purpouses.
#
#Function: TestFiles()
#
#Value: Random testfiles for the pipeline functions including genotypic and map files. First column
#		in genotype file is the names (ids) of samples. Columns 2-6 are arbitary but are there to follow
#		normal ALLELIC genotypic output file structure.
#
#n					The number of samples in the files.
#
#p					The number of markers in the files.
#
#Phenotypes			Vector of names for phenotypes.
#
#MarkerName			Name for markers. The whole name will be MarkerName_n.
#
#chorN				Number of created chromosomes (default is one). 
#
#chromLengths		Vector of chromosome lengths. The total length must be the same as the number of markers. 
#					If there are three chromosomes with markers from 1:200,
#					201:300 and 301:350, the vector should be c(200, 100, 50).
#
#OutDir				Name or location for the output directory.
#
#encoding			"ALLELIC" or "RAW". Plink has 6 columns before genotypic data and allelic encoding
#					RAW has one column before genotypic data (ids), an epty cell in location [1,1] followed by
#					marker names in the first row and the allelic encoding is as dosage (0,1,2).
#					
#BadData			Will the genotypic file contain bad data? Monomorphic alleles or MAF under 0.05 frequency
#					can be created as well as missing data. Vector of length three for c(Monomorphic, MAF, missing
#					values). Zero and values above one mean no and values between mean the frequency. Also, missing
#					values will only be created if encoding option is "ALLELIC" because missing values in "RAW" will
#					make the BGLR function unusable. The missing values will be encoded as "0 0".
#
#
TestFiles<-function(n, p, MarkerName="marker", chromN=1, chromLengths=c(0), OutDir="output", encoding="ALLELIC", BadData=c(0,0,0)){
	dir.create(OutDir, showWarnings=FALSE)
	cat("Output directory created", "\n")
	if(.Platform$OS.type == "windows"){flush.console()}
	setwd(OutDir)
	
	cat("Creating genotypic file", "\n")
	if(.Platform$OS.type == "windows"){flush.console()}
	Gfile<-matrix(data=rep(0) ,nrow=n, ncol=(p+6))
	for(i in 7:(p+6)){
		Gfile[,i]<-sample(x=c(0,1,2), size=n, replace=TRUE)
	}
	Gfile[,1]<-c(100001:(100000+n))
	
	
	if(BadData[1] > 0 & BadData[1] < 1){
		cat("Created monomorphic alleles to random positions", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
		k<-round(p*BadData[1], digits=0)
		for(i in c(sample(c(7:p+6), size=k, replace=FALSE))){
			a<-sum(Gfile[,i]==0, na.rm=TRUE)
			b<-sum(Gfile[,i]==2, na.rm=TRUE)
			if(a>b){
				Gfile[,i]==0
			}
			if(a<b){
				Gfile[,i]==2
			}
			
		}
	}
	
	
	if(BadData[2] > 0 & BadData[2] < 1){
		cat("Created MAF<0.5 to random positions by randomly changing alleles to different values", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
		k<-round(p*BadData[1], digits=0)
		for(i in c(sample(c(7:p+6), size=k, replace=FALSE))){
			a<-sum(Gfile[,i]==0, na.rm=TRUE)
			b<-sum(Gfile[,i]==2, na.rm=TRUE)
			C<-sum(Gfile[,i]==1, na.rm=TRUE)
			maf<-((2*b+C)/(a+b+C))
			for(j in c(sample(c(1:n), size=n, replace=FALSE))){
				if(Gfile[j,i]==2 | Gfile[j,i]==1){
					if(maf>=0.5){
						Gfile[j,i]<-0
						a<-sum(Gfile[,i]==0, na.rm=TRUE)
						b<-sum(Gfile[,i]==2, na.rm=TRUE)
						C<-sum(Gfile[,i]==1, na.rm=TRUE)
						maf<-((2*b+C)/(a+b+C))
					}	
				}
			}
		}
	}
	
		
	if(encoding=="ALLELIC"){
		cat("Changed the encoding to ALLELIC by random nucleotide allocation", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
		for(i in 7:(p+6)){
			alleles=c(0,0)
			alleles<-sample(x=c("A", "T", "G", "C"), size=length(alleles), replace=FALSE)
			zero<-paste(alleles[1], alleles[1])
			one<-paste(alleles[1], alleles[2])
			two<-paste(alleles[2], alleles[2])
			Gfile[Gfile[,i]==0, i]<-zero
			Gfile[Gfile[,i]==1, i]<-one
			Gfile[Gfile[,i]==2, i]<-two
		}
		if(BadData[3] > 0 & BadData[3] < 1){
			cat("Created missing values", "\n")
			if(.Platform$OS.type == "windows"){flush.console()}
			Sum<-round(n*p*BadData[3], digits=0)
			while(Sum > 0){
				a<-sample(x=c(7:(p+6)), size=1)
				b<-sample(x=c(1:n), size=1)
				if(Gfile[b,a]!="0 0"){
					Gfile[b,a]<-"0 0"
					Sum=Sum-1
				}
			}
		}
	}
	
	cat("Making the map file", "\n")
	if(.Platform$OS.type == "windows"){flush.console()}
	map<-matrix(rep(0), nrow=p, ncol=4)
	temp<-0
	for(i in 1:chromN){
		temp<-append(temp, c(rep(i, chromLengths[i])))
	}
	temp<-temp[-1]
	map[,1]<-temp
	temp<-0
	for(i in 1:p){
	temp<-append(temp,paste(MarkerName, "_", i, sep=""))
	}
	temp<-temp[-1]
	map[,2]<-temp
	temp<-0
	for(i in 1:chromN){
	temp<-append(temp, c(1:chromLengths[i]))
	}
	temp<-temp[-1]
	map[,4]<-temp
	
	if(encoding=="RAW"){
		Gfile<-Gfile[,-2:-6]
		temp<-c("	", as.character(map[,2]))
		Gfile<-rbind(temp, Gfile)
	}
	cat("Writing output to", OutDir, "\n")
	if(.Platform$OS.type == "windows"){flush.console()}
	write.table(Gfile, file="genotypes.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
	write.table(map, file="map.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
	cat("All Done!", "\n")
	
	
}
#Example
#TestFiles(n=35, p=350, chromN=3, chromLengths=c(20,10,5), BadData=c(0.0,0.0,0.0), encoding="RAW")

#Step 0.2
#
#Creating simulated phenotypes for testing purpouses.
#
#function: TestPhen()
#
#value: simulated phenotypes with markers in the genotypic model in random positions. 
#
#genotype			Name or location of the genotypic file in RAW -format.
#
#map				Name or location of the map file
#
#Nmarkers			Number of markers in the genotypic model (markers with genotypic effect)
#
#Name				Name for the markers in the genotypic model. The whole names will be 
#					Name_1,...,Name_Nmarkers.
#
#Times				Number of created phenotypes.
#
#VarP				Phenotypic variance
#
#H2					Heritability or heritability range for created phenotypes as vector
#
#uniform			Will the effects in a single phenotype follow a uniform distribution? If not, then
#					Heritability must be a range vector and it is treated as the range from which effects
#					are sampled to each marker. The total actual H2 will then be allowed to vary in each 
#					phenotype.
#
#output				Names of the output phenotypic file, map file, result file (.rda) and the output folder as vector. If
#					fourth value is zero the output is printed in the working directory
#
#

TestPhen<-function(genotype, map, Nmarkers=1, Name="Signf", Times=1, VarP=0, H2=0, uniform=TRUE, output=c("phenotypes","map","results",0)){
	T<-getwd()
	GT<-read.table(genotype, sep="\t", row.names=1, header=TRUE)
	map<-read.table(map, sep="\t", colClasses=c("numeric","character","numeric","numeric"))
	Names<-names(GT)
	PT<-row.names(GT)
	variances<-NULL
	Effects<-NULL
	markernames<-NULL
	Gvar=VarP*H2
	
	
	for(i in 1:Nmarkers){
		variances<-append(variances, paste(Name, "_", i, Sep=""))
		Effects<-append(Effects, paste(Name, "_", i, Sep=""))
		markernames<-append(markernames, paste(Name, "_", i, Sep=""))
	}
	for(i in 1:Times){
		markerN<-sample(x=c(1:ncol(GT)), size=Nmarkers)
		if(length(H2)>1){
			if(uniform==TRUE){
				effectVar<-rep(round(runif(1, min=Gvar[1], max=Gvar[2]), digits=4), Nmarkers)
				} else {
				effectVar<-round(runif(Nmarkers, min=Gvar[1], max=Gvar[2]), digits=4)
				}
			
		} else {
			effectVar<-rep(Gvar, Nmarkers)
		}
		if(uniform==TRUE){
		effectVar<-effectVar/Nmarkers
		}
		fenot<-rnorm(n=nrow(GT), mean=0, sd=sqrt(VarP)) 
		variances<-cbind(variances, effectVar)
		colnames(variances)[(i+1)]<-paste("Phenotype", i, sep="")
		
		for(k in 1:length(markerN)){
			if(grepl(pattern=Name, map[markerN[k],2])==TRUE){
				map[markerN[k],2]<-paste(map[markerN[k],2], "&", markernames[k], sep="")
				} else {
				map[markerN[k],2]<-markernames[k]
			}
		}
		effectT<-rep(0,Nmarkers)
		for(j in 1:length(markerN)){
	
			Pbi<-sum(GT[,markerN[j]])/nrow(GT)
			Phet<-length(which(GT[,markerN[j]]==1))/nrow(GT)
			Phom<-length(which(GT[,markerN[j]]==2))/nrow(GT)
	
			effect<-sqrt(effectVar[j]/(Phet-2*Phet*Pbi+4*Phom-4*Phom*Pbi+Pbi^2)) 
			fenot[GT[,markerN[j]]==2] <- fenot[GT[,markerN[j]]==2] + (effect*2)
			fenot[GT[,markerN[j]]==1] <- fenot[GT[,markerN[j]]==1] + (effect)
			names(GT)[markerN[j]]<-markernames[j]
			effectT[j]<-effect
		}
		Effects<-cbind(Effects ,effectT)
		colnames(Effects)[(i+1)]<-paste("Phenotype", i, sep="")
		PT<-cbind(PT, fenot)
		colnames(PT)[(i+1)]<-paste("Phenotype", i, sep="")
	}
	if(output[4]!=0){
		
		dir.create(file.path(output[4]), showWarnings = FALSE)
		setwd(file.path(Output[4]))
	
	}
	
	Results<-list(Phenotypes=PT, Effects=Effects, Variances=variances, map=map)
	save(Results, file=paste(output[3], ".rda", sep=""))
	write.table(PT, file=output[1], sep="\t", quote=FALSE, row.names=FALSE)
	write.table(map, file=output[2], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	setwd(T)
}

#Example
#TestPhen(genotype="genotypes.txt", map="map.txt", Nmarkers=4, Name="Signf", Times=10, VarP=1, H2=c(0.3, 0.4),uniform=FALSE, output=c("phenotypes","map","results",0))



#Step 1
#Removing unwanted chromosomes from the data and quality controll. This function will remove unwanted chromosomes
#and check the datasets for missing information. Markers and samples are checked for bad call rate. Minor allele
#frequency can be controlled as well. Markers in same genomic positions can be cleaned. Tests fo heterozygozity
#and for Hardy-Weinberg equilibrium can be made as well. The output print will include information about the 
#removed data. The controlled datasets cand be separated into chromosomes or saved into output directory as 
#the whole genome or both options can be used.
#
#function: qualityCtrl()
#
#Value: Quality controlled dataset(s). A print containing information about the removed data and tests made.
#
#genotypes			Name of the genotype file (tab delimited)
#
#map				Name of the map file (tab delimited)
#
#Filter				Vector of chromosomes to be removed or NULL
#
#Prune				Will the data be checked for bad call rate and those 
#					positions removed. TRUE or FALSE.
#
#sampleCallRate		Will the samples (animals) be checked for bad call rate
#					and those samples removed? Value of 0 or higher than 1 mean no, values above 0
#					(but smaller than 1) are considered as the screening limit. 
#
#MAF				Will the markers be checked for minor allele frequency (MAF) and will those markers
#					be removed? Value of 0 or higher than 1 mean no, values above 0
#					(but smaller than 1) are considered as the screening limit.
#
#CleanDouble		Will the markers be checked for duplicate genomic positions? If TRUE then the other markers 
#					position will increased by one base.
#
#Tests				Will test for heterozygozity and H-W equilibrium be made? 
#
#Separate			If TRUE the datasets will be returned separated into chromosomes. If FALSE datasets 
#					will be returend as genomic. Value can also be BOTH for both options to be used.
#
#Output				A vector of three values. First is the name for the output chromosome file(s), the second 
#					is a name for the output map file(s), and the third is the name (or location) for the output 
#					directory.
#

qualityCtrl<-function(genotypes, map, Filter=NULL,  prune=TRUE, sampleCallRate=0, MAF=0.05, cleanDouble=TRUE, Tests=TRUE, separate=TRUE, output=c("chromosome", "map", "output")){
	
	markers<-read.table(genotypes, sep="\t", colClasses="character", check.names=FALSE)
	map<-read.table(map, sep="\t", colClasses="character", check.names=FALSE)
	cat("data tables read", "\n")
	cat("Starting with", ncol(markers[,-1:-6]), "markers in ", nrow(markers), "samples.", "\n")
	removedData<-list(pruned=0, removedSamples=0, MAF_removed=0, duplicates=0)
	results<-list(removedData=0, tests=0)
	if(length(output)>2){
	dir.create(output[3], showWarnings=FALSE)
	cat("Output directory created", "\n")
	setwd(output[3])
	}
	if(.Platform$OS.type == "windows"){flush.console()}
	if(is.null(Filter)==TRUE){
		temp<-list(markers,map)
		if(length(temp)==2){
			cat("Data read into a list", "\n")
		}
	} else {
		if(is.vector(Filter)==FALSE | is.numeric(Filter)==FALSE){
			cat("Error: Filter must be either NULL or numeric vector.", "\n")
			break
		}
				
		mapc<-c("d1","d2","d3","d4","d5","d6", c(map[,1]))
		markers<-markers[,!(mapc %in% Filter)]
		map<-map[!(map[,1] %in% Filter),]
		temp<-list(markers,map)
		if(length(temp)==2){
			cat("Data read into a list and chromosomes ", Filter, " filtered", "\n")
		}
	}
	
		
	if(prune==TRUE){
		cat("Pruning data...", "\n")
		a=0.05
		b=nrow(temp[[1]])
		c=a*b
		index=0
		for(i in 7:ncol(temp[[1]])){
			missing<-which(temp[[1]][,i]=="0 0")
			if(length(missing)>c){
				index<-c(index,i)
			}
		}
		
		if(length(index)>1){
			index<-index[-1]
			temp[[1]]=temp[[1]][,-index]
			index<-(index-6)
			removedData$pruned<-temp[[2]][index,]
			temp[[2]]=temp[[2]][-index,]
		}
		n<-(ncol(temp[[1]])-6)
		
		cat("Markers checked for bad call rate.", length(index), "markers removed", "\n")
		cat(n, "markers left", "\n")
	}	
	if(sampleCallRate>0 & sampleCallRate<1){
		cat("Cheking samples for bad call rate", "\n")
		index<-0
		n<-ncol(temp[[1]][,-1:-6])
		a<-c(seq(1,n,by=5000))
		b<-c(seq(0,n,by=5000), n)
		b<-b[-1]
		d<-1
		for(j in 1:(length(a))){
			
			dat<-temp[[1]][,a[j]:b[j]]
			Gz<-matrix(data=rep(0), nrow=nrow(dat), ncol=ncol(dat))
			for(i in 1:ncol(dat)){
				Gz[,i][dat[,i]=="0 0"]<-1
			}
			if(j==1){
				sums<-rowSums(Gz)
			} else {
				sums<-sums+rowSums(Gz)
			}
			cat(b[j], "markers checked", "\n")
			if(.Platform$OS.type == "windows"){flush.console()}
		}
		cat("Calculating missing marker sums for samples...", "\n")
		sums<-1-(sums/(ncol(markers)))
		removedData$removedSamples<-temp[[1]][sums<=sampleCallRate,1]
		temp[[1]]<-temp[[1]][sums>sampleCallRate,]
		index<-length(removedData$removedSamples)
		cat("Sample call rate checked", index, "samples removed", "\n")
		n<-nrow(temp[[1]])
		cat(n, "samples left", "\n")
	}

	if(MAF>0 & MAF<1){
		cat("MAF check starting...", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
		index<-0
		d<-0
		a<-0
		mono<-0
		maf<-0
		n<-nrow(temp[[1]])
		Gz1<-matrix(data=rep(0), nrow=5, ncol=n)
		for(i in 7:ncol(temp[[1]])){
			Gz<-Gz1
			Gz[5,][temp[[1]][,i]!="A A" & temp[[1]][,i]!="C C" & temp[[1]][,i]!="G G" & temp[[1]][,i]!="T T"]<-1
			Gz[5,][temp[[1]][,i]=="0 0"]<-0
			Gz[1,][temp[[1]][,i]=="A A"]<-2 ; Gz[2,][temp[[1]][,i]=="C C"]<-2 ; Gz[3,][temp[[1]][,i]=="G G"]<-2 ; Gz[4,][temp[[1]][,i]=="T T"]<-2
			vec<-rowSums(Gz)
			minor<-((min(vec[1:4][vec[1:4]!=0])+vec[5])/sum(vec))
			major<-((max(vec[1:4])+vec[5])/sum(vec))
			
			if(minor==major){
				if(length(vec[1:4][vec[1:4]!=0])>1){
					d<-c(d,i)
					next
					} else {
						if(sum(vec[5])==0){
							mono<-c(mono,i)
							index<-c(index,i)
						} else {
							index<-c(index,i)
							a<-c(a,i)
						}
					}
				} else {
					if(minor<MAF){
						index<-c(index,i)
						maf<-c(maf,i)
					}
				}
			
			if((i %% 5000)==0 | i==ncol(temp[[1]])){
				cat(i, "markers checked", "\n")
				if(.Platform$OS.type == "windows"){flush.console()}
			}
		}
	
		if(length(index)>1){
			index<-index[-1]
			if(length(d)>1){d<-(length(d)-1)}
			if(length(a)>1){a<-(length(a)-1)}
			if(length(mono)>1){mono<-(length(mono)-1)}
			if(length(maf)>1){maf<-(length(maf)-1)}
			temp[[1]]<-temp[[1]][,-index]
			index<-(index-6)
			removedData$MAF_removed<-temp[[2]][index,]
			temp[[2]]=temp[[2]][-index,]
			index<-length(index)
		}
		cat("MAF checked. ", maf, "markers removed.", "\n")
		cat(mono, "markers removed for monomorphicity.", "\n", a ,"markers removed for having no minor allele homozygotes.", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
	}
	
	if(cleanDouble==TRUE){
		cat("cleaning map from duplicates. (duplicated position is now +(ni*1bp)", "\n")
		dat<-temp[[2]][,4][duplicated(temp[[2]][,4])]
		v<-c(1:4)
		for(i in 1:length(dat)){
			fm<-temp[[2]][temp[[2]][,4] %in% dat[i],]
			fm1<-rownames(fm)
			fm2<-as.numeric(fm[,4])
			for(j in 0:(length(fm1)-1)){
				fm2[j+1]<-fm2[j+1]+j
			}
			temp[[2]][fm1,4]<-fm2
			v<-rbind(v,fm)
		}
		removedData$duplicates<-v
	}
	
	results[[1]]<-removedData
	tests<-0
	results[[2]]<-tests
	cat(ncol(temp[[1]][,-1:-6]), "markers left", "\n")
	
	if(Tests==TRUE){		
		cat("Heterozygozity test for markers starting.", "\n")
		cat("Hardy-Weinberg test for markers starting.", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
		gen<-c(rep.int(0, nrow(temp[[2]])))
		tests<-list(markerMap=temp[[2]],n=gen, P=gen, Q=gen, Obs=data.frame(gen,gen,gen), Exp=data.frame(gen,gen,gen),  HetObs=gen, HetExp=gen )
		names(tests$Obs)<-c("PP", "PQ", "QQ")
		names(tests$Exp)<-c("PP", "PQ", "QQ")
		n<-nrow(temp[[1]])
		Gz1<-matrix(data=rep(0), nrow=5, ncol=n)
		for(i in 7:ncol(temp[[1]])){
			j=i-6
			Gz<-Gz1
			Gz[5,][temp[[1]][,i]!="A A" & temp[[1]][,i]!="C C" & temp[[1]][,i]!="G G" & temp[[1]][,i]!="T T"]<-1
			Gz[5,][temp[[1]][,i]=="0 0"]<-0
			Gz[1,][temp[[1]][,i]=="A A"]<-2 ; Gz[2,][temp[[1]][,i]=="C C"]<-2 ; Gz[3,][temp[[1]][,i]=="G G"]<-2 ; Gz[4,][temp[[1]][,i]=="T T"]<-2
			vec<-rowSums(Gz)
			if(length(vec[vec!=0])==3){										
				tests$P[j]<-Pi<-((max(vec[1:4])+sum(vec[5]))/sum(vec))				#Pi
				tests$Q[j]<-Qi<-((min(vec[1:4][vec[1:4]!=0])+sum(vec[5]))/sum(vec)) #Qi
				tests$Obs[j,1]<-max(vec[1:4])										#ObsPPi
				tests$Obs[j,3]<-min(vec[1:4][vec[1:4]!=0])							#ObsQQi
				tests$Obs[j,2]<-sum(vec[5])											#ObsQPi
				tests$n[j]<-N<-sum(colSums(Gz!=0))									#n
				tests$Exp[j,3]<-N*Qi^2												#ExpQQi
				tests$Exp[j,1]<-N*Pi^2												#ExpPPi
				tests$Exp[j,2]<-2*Pi*Qi*N											#ExpQPi
				tests$HetObs[j]<-sum(vec[5])/N										#HetObs
				tests$HetExp[j]<-2*Pi*Qi											#HetExp
			}
			if(length(vec[vec!=0])<3){
				tests$P[j]<-Pi<-((max(vec[1:4])+sum(vec[5]))/sum(vec))
				tests$Q[j]<-Qi<-sum(vec[5])
				tests$Obs[j,1]<-max(vec[1:4])
				tests$Obs[j,2]<-sum(vec[5])
				tests$n[j]<-N<-sum(colSums(Gz!=0))
				tests$Exp[j,3]<-N*Qi^2
				tests$Exp[j,1]<-N*Pi^2
				tests$Exp[j,2]<-2*Pi*Qi*N
				tests$HetObs[j]<-sum(vec[5])/N
				tests$HetExp[j]<-2*Pi*Qi
			}
			if((i %% 5000)==0 | i==ncol(temp[[1]])){
				cat(i, "markers checked", "\n")
				if(.Platform$OS.type == "windows"){flush.console()}
			}
		}
		results[[2]]<-tests
	}		
		
	cat("All tests completed. Writing data", "\n")	
	
			
	if(separate==TRUE | separate=="BOTH"){
		
		mapc<-c("d1","d2","d3","d4","d5","d6", c(temp[[2]][,1]))
		chromnumbers<-unique(temp[[2]][,1])
		idcolumns<-temp[[1]][,1:6]
		for(i in 1:length(chromnumbers)){
			mapc[1:6]<-chromnumbers[i]
			name<-paste(output[1],chromnumbers[i], sep="")
			nameMap<-paste(output[2], chromnumbers[i], sep="")
			write.table(temp[[1]][,mapc %in% chromnumbers[i]], file=name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
			write.table(temp[[2]][temp[[2]][,1] %in% chromnumbers[i],], file=nameMap,  sep="\t",col.names=FALSE, row.names=FALSE, quote=FALSE)
		}
		
		cat ("Data separated into chromosomes and written in", getwd(), "\n")
		if(separate!="BOTH"){
			setwd("..")
		}
				
			
	}	
	
	if(separate==FALSE | separate=="BOTH"){
			
			write.table(temp[[1]], file=output[1], sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
			write.table(temp[[2]], file=output[2], sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
			cat("Data written into files in", getwd(), "\n")
			setwd("..")
				
			
			
			
	}
	return(results)	
}		
		
#Example:
#data<-qualityCtrl("markers.ped", "map.map", Filter=c(0,1,2),  prune=TRUE, sampleCallRate=0.90, MAF=0.05, cleanDouble=TRUE, Tests=TRUE, separate="BOTH", output=c("chromosome", "map", "output"))

		
#Step 3.
#Imputing and phasing the data using SHAPEIT

#See SHAPEIT manual for further details if needed. (www.shapeit.fr/)
# 

#Example of SHAPEIT commandline in PERL below:
#shapeit --input-ped chromosomefile mapfile \
#--output-max chromosome.haps chromosome.sample

#CSC taito shell has module for SHAPEIT which can be launched
#by command:
#module load shapeit
#
#Because the phasing and imputing can be very time consuming, in CSC environment 
#a batch job should be launched for example by using the batch job wizard.
#The phasing might not run to completion if not.


#Step 4. 
#Changing the SHAPEIT phased data into "RAW" format. This function will scan the working directory 
#for files ending in "N.haps", where the N is the chromosome number and then transform them into "RAW" 
#format. The original files will not be destroyed. 
#
#function: PhasedToBack() 
#
#Value: genotypic files in "RAW" (.BGLR) or "ALLELIC" (.ped) format.
#
#chrnumbs				The chromosome files that are transformed. Argument is a vector of
#						chromosome numbers to be changed.
#
#samples				A vector of sample names (Ids)
#
#


PhasedToBack<-function(chrnumbs=NULL, samplesOrIdS=0, ALLELIC=TRUE){
	if(is.null(chrnumbs)==TRUE){
	break
	}
	files<-list.files(pattern=".haps", full.names=TRUE)
	files2<-list.files(pattern=".haps")
	
	
	repeat {
		
		filename<-files[grep(paste(tail(chrnumbs, n=1), ".haps", sep=""), files)]
		filename2<-files2[grep(paste(tail(chrnumbs, n=1), ".haps", sep=""), files2)]
		if(ALLELIC==TRUE){
			writtenin<-gsub(".haps", ".ped", filename2, fixed=TRUE)
		}
		if(ALLELIC==FALSE){
			writtenin<-gsub(".haps", ".BGLR", filename2, fixed=TRUE)
		}
		
		haplot<-read.table(filename, sep=" ")
		markerIDs<-c("	", as.character(haplot[,2]))
		newdata<-matrix(rep(0), nrow=(ncol(haplot[,-1:-5])/2), ncol=nrow(haplot))
		if(ALLELIC==TRUE){
			alleles<-haplot[,4:5]
		}
		haplot<-haplot[,-1:-5]
		haplot<-as.matrix(haplot)
		
		cat("Chromosome", tail(chrnumbs,1), "next", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
		
		
		a<-seq(1,ncol(haplot),2)
		b<-a+1
		s<-nrow(newdata)
	
		if(tail(b,1)!=ncol(haplot)){
			cat("Something has failed horribly")
			}
					
		for(l in 1:length(a)){
			vec<-rowSums(haplot[,a[l]:b[l]])
			newdata[l,]<-vec
		}
		
		if(ALLELIC==TRUE){
			for(i in 1:ncol(newdata)){
				newdata[,i][newdata[,i]==0]<-paste(alleles[i,1], alleles[i,1], sep=" ")
				newdata[,i][newdata[,i]==1]<-paste(alleles[i,1], alleles[i,2], sep=" ")
				newdata[,i][newdata[,i]==2]<-paste(alleles[i,2], alleles[i,2], sep=" ")
			}
		}
		
		if(ALLELIC==FALSE){
			for(i in 1:ncol(newdata)){
				if(sum(newdata[,i])>s){
					newdata[,i][newdata[,i]==2]<-3
					newdata[,i][newdata[,i]==0]<-2
					newdata[,i][newdata[,i]==3]<-0
				}
			}
		}
		if(ALLELIC==TRUE){
			newdata<-cbind(samplesOrIdS, newdata)
		}
				
		if(ALLELIC==FALSE){
			newdata<-cbind(samplesOrIdS, newdata)
			newdata<-rbind(markerIDs, newdata)
		}
		
		write.table(newdata, file=writtenin, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
		
		cat("Chromosome", tail(chrnumbs,1), "transformed", "\n")
		if(.Platform$OS.type == "windows"){flush.console()}
		
		chrnumbs<-chrnumbs[-tail(chrnumbs, n=1)]
		files<-files[!files %in% filename]
		files2<-files2[!files2 %in% filename2]
		
		if(length(chrnumbs)==0){
			break
		}
	}	
}

 
#Example:
#sample vector can be extracted for example from the first column of the file containing the ALLELIC genotypes.
#PhasedToRaw(chrnumbs=c(1:3), samples=samples, ALLELIC=TRUE)
 
 
#Step 5.
#Single marker analysis for the SIS prosedure. For SIS prosedure the single marker method -analysis was 
#done by the program TASSEL (method GLM) in this pipeline, but it could be done with
#any single marker method. Genotypic and Map files can be loaded into tassel in the same
#ALLELIC-format used in this pipeline (choose "load plink" from the "data/load/choose filetype to load" options). 
#For the phenotype file, two more lines needs to be added. A line containing the indicator: <Phenotype> and another
#Line containing the indicators for taxa (id) and data columns. You can do this manually or use the following code:
#
#Function: AddLines()
#
#Value: Adds the needed indicator lines to the phenotype file when doing single marker analysis with Tassel
#		and returns the file as an object.
#
#
#Phenotypes					The phenotypic file as an object
#
#

AddLines<-function(Phenotypes){
	n=ncol(Phenotypes)
	p=append(paste("<Phenotype>"), rep(paste("	", sep="\t"), (n-1)))
	p=rbind(p,append(paste("taxa"), rep(paste("data", sep="\t"), (n-1))))
	Phenotypes=rbind(p, Phenotypes)
	return(Phenotypes)
}

#Step 6.
#Selection for the results of single marker analysis
#
#function: PreSelection()
#
#Value: pre selected markers as a new dataset (genotype and map) printed into working directory
#
#n					Number of markers selected
#
#genotype			Name or location for the Genotypic file in "RAW"-format
#
#map				Name or location for the Map file
#
#TR					Name or location for the Tassel results file. This can contain more than
#					one phenotype, they will be scanned and separated automatically.
#

PreSelection<-function(n, genotype, map, TR){
	if(is.numeric(n)==FALSE | n%%1!=0){cat("n must be numeric integer"); break}
	testfile<-read.table(TR, sep="\t", dec=",", header=TRUE, colClasses=c(rep("character", 2), rep("numeric", 15)))
	G<-read.table(genotype, sep="\t", row.names=1, header=TRUE)
	map<-read.table(map, sep="\t", colClasses=c("numeric","character","numeric","numeric")))
	types<-table(testfile[,1])
	Ntypes<-length(types)
	for(i in 1:Ntypes){
		s<-names(types[i])
		Data<-testfile[testfile[,1]==s,]
		p<-Data$p
		p<-order(p, decreasing=FALSE)
		p<-p[1:n]
		p<-sort(b)
		temp<-G[,b]
		Tmap<-map[b,]
		names(temp)<-Tmap[,2]
		
		write.table(temp, file=paste(names(types[i], n, "_genotype", sep="")) sep="\t", quote=FALSE)
		write.table(Tmap, file=paste(names(types[i], n, "_map", sep="")) sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
	
}

#Step 7
#The modified BGLR function. This is just the modiefied version of the BGLR function for the BGLR package.
#Load the BGLR package into R first and then replace the original function with this one (by either copying
#and pasting from this file or loading from a separate file). The getResults() function is the main part of
#the code added inside the BGLR() function.
#
#function: BGLR(), getResults()
#
#Added arguments:
#
#BigMatrixName				The name for indicatormatrix to be stored in results. If NULL, the indicator matrix will not be stored
#
#map						The map file as an object
#
#Windowed					Will a windowed analysis also be performed. The windows will be very short, or only the adjacent markers.
#
#
#
#
#Additional Notes for modified version:
#Modified version of BGLR and its subfunctions. Modifications are in the setLT.BayesBandC function and BGLR function for BayesB and BayesC model prints.
#Usage: load BGLR package with library(BGLR) and then this modified script with source("ThisFile.r"). This script will also
#Check for package bigmemory and install/load it as needed. There is also one more argument added to BGLR function, BigMatrixName
#this is used to define the name/path for saving the indicator matrix.
#Other new outputs in modified version:
#$ETA$MRK$BF 		= 	Posterior Predictive Bayes factors for covariates in genomic order
#$ETA$MRK$MCMCinD	=	Model selection indicators summed over iterations
#$ETA$MRK$Pi		=	Values of  π parameter (priors) over iterations as vector


getResults<-function(fm=fm, Windowed=TRUE, indicatorMatrx=NULL, map=map, BFlimit=SignfLim, BFlimit2=SignfLimKa, BFlimitCustom=CustomLimit, BF=BF, b=b, Pi=Pi) {
	
	windowed<-list(overL=0, overT=0, BFlimit=0, overTka=0, BFlimitMean=0, overTcustom=0, BFlimitCustom=BFlimitCustom)
	results<-list(overL=0, overT=0, BFlimit=0, overTka=0, BFlimitMean=0, overTcustom=0, BFlimitCustom=BFlimitCustom)
	Xij<-list(results=0, windowed=0, BFw=0)
	Names<-names(b)
	Names<-gsub(".", "-", Names, fixed=TRUE)
	overL<-map[(map[,2] %in% Names),]
	if(is.null(indicatorMatrx)==FALSE){
		Mx<-indicatorMatrx
	}
	OverL<-cbind(overL, BF)
	names(BF)<-Names
	#95% two-way test for BF components
	
		getmode <- function(v) {
			uniqv <- unique(v)
			uniqv[which.max(tabulate(match(v, uniqv)))]
		}
		BFmod<-0
		
		inc1<-0
		inc2<-0
		inc3<-0
		for(i in 1:ncol(Mx)){
			M<-nrow(Mx)
			pI1<-(sum(Mx[,i])/M)
			pI1O=pI1/(1-pI1)
			Ipi<-Mx[,i]*(Pi/(1-Pi))
			Ipi<-Ipi[Ipi>0]
			Ipi<-pI1O/Ipi
			fr1<-length(which(Ipi>BFlimit))/length(Ipi)
			fr2<-length(which(Ipi>BFlimit2))/length(Ipi)
			fr3<-length(which(Ipi>BFlimitCustom))/length(Ipi)
			BFmod<-append(BFmod, getmode(Ipi))
						
			inc1<-append(inc1, fr1)
			inc2<-append(inc2, fr2)
			inc3<-append(inc3, fr3)
		}
		
		Test<-inc1[-1]
		Test2<-inc2[-1]
		Test3<-inc3[-1]
		
		BFmod<-BFmod[-1]
		overT<-BF[which(Test>0.95)]
		overTka<-BF[which(Test2>0.95)]
		overTcustom<-BF[which(Test3>0.95)]
		OverL<-cbind(OverL, BFmod, Test)
	
	
	OverL<-OverL[OverL$BF>BFlimit2,]
	names(OverL)<-c("Chrom", "Marker", "zero", "bp", "BF", "BFmode", "Test")
	results$overL<-OverL
	results$overT<-overT
	results$overTka<-overTka
	results$overTcustom<-overTcustom
	results$BFlimit<-BFlimit
	results$BFlimitMean<-BFlimit2
	
	if(Windowed==TRUE){
		
		M<-nrow(Mx)
		BFmod<-0
		BFw=0
		inc1<-0
		inc2<-0
		inc3<-0
		BFw<-0
		for(i in 1:ncol(Mx)){
			d<-0
			if(i!=1 & i!=ncol(Mx)){
				d = Mx[,i] + Mx[,(i-1)] + Mx[,(i+1)]
				d[d>1]<-1
			}	
			if(i==1){
				d = Mx[,i] + Mx[,(i+1)]
				d[d>1]<-1
			}
			if(i==ncol(Mx)){
				d = Mx[,i] + Mx[,(i-1)]
				d[d>1]<-1
			}
			pI1<-(sum(d)/M)
			pI1O=pI1/(1-pI1)
			Ipi<-d*(Pi/(1-Pi))
			Ipi<-Ipi[Ipi>0]
			Ipi<-pI1O/Ipi
			fr1<-length(which(Ipi>BFlimit))/length(Ipi)
			fr2<-length(which(Ipi>BFlimit2))/length(Ipi)
			fr3<-length(which(Ipi>BFlimitCustom))/length(Ipi)
			BFmod<-append(BFmod, getmode(Ipi))
			BFw<-append(BFw, mean(Ipi))
			
			inc1<-append(inc1, fr1)
			inc2<-append(inc2, fr2)
			inc3<-append(inc3, fr3)
			
		
		}
		Test<-inc1[-1]
		Test2<-inc2[-1]
		Test3<-inc3[-1]
		
		BFmod<-BFmod[-1]
		BFw<-BFw[-1]
		OverL<-cbind(overL, BFw)
		names(BFw)<-Names
		overT<-BFw[which(Test>0.95)]
		overTka<-BFw[which(Test2>0.95)]
		overTcustom<-BFw[which(Test3>0.95)]
		overL<-cbind(overL, BFw)
		overL<-cbind(overL, BFmod, Test)
		
		#test for BF limit
		overL<-overL[overL$BFw>BFlimit2,]
		names(overL)<-c("Chrom", "Marker", "zero", "bp", "BFw", "BFmode", "Test")
		windowed$overL<-overL
		windowed$overT<-overT
		windowed$overTka<-overTka
		windowed$overTcustom<-overTcustom
		windowed$BFlimit<-BFlimit
		windowed$BFlimitMean<-BFlimit2
		
	
	
	}
	
	Xij$results<-results
	Xij$windowed<-windowed
	Xij$BFw<-BFw
	return(Xij)
	
}
#End of notes for modified version.




#This function creates an incidence matrix that will be included in the 
#linear term of the model
#Arguments: LT, Linear term, an object of the class "formula" that also includes
#optionally a data.frame to obtain the information
#It returns the incidence matrix
set.X=function(LT)
{	
	flag=TRUE
	n_elements=length(LT)
	i=0
	while(i<=n_elements & flag)
	{
	   i=i+1;
	   if(class(LT[[i]])=="formula")
	   {
	   		flag=FALSE
			rhs=LT[[i]]
			if(is.null(LT$data))
			{
				mf = model.frame(formula=rhs)
			}else{
				mf = model.frame(formula=rhs,data=LT$data)
			}
    		X = model.matrix(attr(mf, "terms"), data=mf)
    		Xint = match("(Intercept)", colnames(X), nomatch=0L)
    		if(Xint > 0L) X = X[, -Xint, drop=FALSE]
	   }
	}
	if(flag) stop("Unable to build incidence matrix, wrong formula or data\n")
	return(X)
}

## Fixed Effects ##################################################################
#Function for initializing regression coefficients for Fixed effects.
#All the arguments are defined in the function BGLR
setLT.Fixed=function(LT,n,j,y,weights,nLT,saveAt,rmExistingFiles,groups,nGroups)
{

    if(is.null(LT$X)) LT$X=set.X(LT)

    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
	
    if(any(is.na(LT$X)))
    { 
	stop(paste(" LP ",j," has NAs in X",sep=""))
    }

    if(nrow(LT$X)!=n)
    {
        stop(paste(" Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))
    }

    #weight inputs if necessary
        
    LT$X=sweep(LT$X,1L,weights,FUN="*")        #weights

    if(!is.null(groups))
    {
        x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
        for(g in 1:nGroups)
        {
                x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
        }
        LT$x2=x2;
    }else{
        LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    }	
    
    #Objects for saving posterior means from MCMC
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)   

    fname=paste(saveAt,LT$Name,"_b.dat",sep="")

    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$fileOut=file(description=fname,open="w")
    tmp=LT$colNames
    write(tmp, ncolumns = LT$p, file = LT$fileOut, append = TRUE)
    LT$X=as.vector(LT$X)
    LT$x2=as.vector(LT$x2)
    LT$varB=1e10
    return(LT)
}

## Gaussian Regression ############################################################
#Function for initializing regression coefficients for Ridge Regression.
#All the arguments are defined in the function BGLR
setLT.BRR=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,groups,nGroups,verbose)
{
    #Check inputs

    if(is.null(LT$X)) LT$X=set.X(LT)
       
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
	
    if(any(is.na(LT$X)))
    { 
      stop(paste(" LP ",j," has NAs in X",sep=""))
    }
    
    if(nrow(LT$X)!=n)
    {
      stop(paste(" Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))
    }   

    #Weight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights

    if(!is.null(groups))
    {
    
	x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
	for(g in 1:nGroups)
	{
		x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
	}
        LT$x2=x2;
    }else{
	LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    }  

    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)

    #Default df for the prior assigned to the variance of marker effects
    if(is.null(LT$df0))
    {
	LT$df0=5

	if(verbose)
	{
		cat(paste(" Degree of freedom of LP ",j,"  set to default value (",LT$df0,").\n",sep=""))
	}
    }

    if(is.null(LT$R2))
    { 
        LT$R2=R2/nLT
    }

 
    #Default scale parameter for the prior assigned to the variance of marker effects
    if(is.null(LT$S0))
    {
        if(LT$df0<=0) stop("df0>0 in BRR in order to set S0\n")

	LT$MSx=sum(LT$x2)/n-sumMeanXSq       
	LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(LT$MSx))*(LT$df0+2)  
	
	if(verbose)
	{
		cat(paste(" Scale parameter of LP ",j,"  set to default value (",LT$S0,") .\n",sep=""))
	}
    }

    
    #Objects for saving posterior means from MCMC
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$varB=LT$S0/(LT$df0+2)
    LT$post_varB=0                 
    LT$post_varB2=0
    
    fname=paste(saveAt,LT$Name,"_varB.dat",sep=""); 
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")
    LT$X=as.vector(LT$X)
    LT$x2=as.vector(LT$x2)

    return(LT)
}

#Ridge regression using sliding windows
#This is just a Ridge Regression with sliding windows, 
#LT has two extra attributes: window_list, and n_windows
#If n_windows is given the the program will obtain a list with the markers in each sliding window.

setLT.BRR_windows=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,verbose)
{

    #Check the inputs
    if(is.null(LT$X)) LT$X=set.X(LT)
   
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)

    if(is.null(LT$windows_list) & is.null(LT$nwindows)) stop("Provide windows_list or nwindows\n");
    if((!is.null(LT$windows_list)) & (!is.null(LT$nwindows))) stop("Provide only windows_list or nwindows but no both\n");
	
    if(any(is.na(LT$X)))
    { 
      stop(paste(" LP ",j," has NAs in X",sep=""))
    }
    
    if(nrow(LT$X)!=n)
    {
      stop(paste(" Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))
    }   

    #Weight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
    

    if(is.null(LT$df0))
    {
	LT$df0=5
	if(verbose)
	{
		cat(paste(" Degree of freedom of LP ",j,"  set to default value (",LT$df0,").\n",sep=""))
	}
    }

    if(is.null(LT$R2))
    { 
        LT$R2=R2/nLT
    }

    if(is.null(LT$S0))
    {
        if(LT$df0<=0) stop("df0>0 in BRR in order to set S0\n")

	LT$MSx=sum(LT$x2)/n-sumMeanXSq       
	LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(LT$MSx))*(LT$df0+2) 
	if(verbose)
	{ 
		cat(paste(" Scale parameter of LP ",j,"  set to default value (",LT$S0,") .\n",sep=""))
	}
    }

    if(is.null(LT$windows_list))
    {
    	windows_list=list()
        nwindows=LT$nwindows

    	if(nwindows>1)
    	{
        	s=as.integer(LT$p/nwindows)
        	for(i in 1:(nwindows-1))
        	{
           		windows_list[[i]]=c(1:s)+s*(i-1)
        	}
        	windows_list[[nwindows]]=c(((nwindows-1)*s+1):LT$p)
    	}else{
          	stop(paste("It does not make any sense to call this function with ",nwindows, " window(s)!!!\n"))
    	}

    	LT$windows_list=windows_list
    }
    
    if(is.null(LT$nwindows))
    {
	LT$nwindows=length(LT$windows_list)
        if(LT$nwindows<=1) stop(paste("The length of the window_list that you provided is ",LT$nwindows," we are expecting a list of length at least 2\n"));
    }
	
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
    LT$post_varB=0                 
    LT$post_varB2=0

    fname=paste(saveAt,LT$Name,"_varB.dat",sep=""); 
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")
    LT$X=as.vector(LT$X)

    return(LT)
}

## Bayesian LASSO ############################################################
## The well known Bayesian LASSO (Park and Casella, 2008) and 
## de los Campos et al (2009). 
#  This functions simply sets hyper-parameters for quantities involved in BL regression

setLT.BL=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,verbose)
{
    #Check the inputs
    if(is.null(LT$minAbsBeta)) LT$minAbsBeta=1e-9

    if(is.null(LT$X)) LT$X=set.X(LT)

    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
		
    if(any(is.na(LT$X)))
    {
       stop(paste("LP ",j," has NAs in X",sep=""))
    }

    if(nrow(LT$X)!=n)
    {
        stop(paste(" Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))
    } 

    #Wheight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)	

    LT$MSx=sum(LT$x2)/n-sumMeanXSq

    # Prior
    if(is.null(LT$R2))
    { 
       LT$R2=R2/nLT
    }
	
    # Setting default value of lambda
    if(!is.null(LT$lambda))
    { 
    	if(LT$lambda<0)
        {
		stop(" lambda should be positive\n")
	}  
    }
    if(is.null(LT$lambda))
    {
		LT$lambda2=2*(1-R2)/(LT$R2)*LT$MSx
		LT$lambda=sqrt(LT$lambda2)
		if(verbose)
		{
			cat(paste(" Initial value of lambda in LP ",j," was set to default value (",LT$lambda,")\n",sep=""))
		}
    }else{
	if(LT$lambda<0) stop(" lambda should be positive\n");
        LT$lambda2=LT$lambda^2
    }
	
    # Checking lambda-type
    if(is.null(LT$type))
    {
		LT$type="gamma"
		if(verbose)
		{
			cat(paste("  By default, the prior density of lambda^2 in the LP ",j,"  was set to gamma.\n",sep=""))
		}
    }else{
		if(!LT$type%in%c("gamma","beta","FIXED")) stop(" The prior for lambda^2 should be gamma, beta or a point of mass (i.e., fixed lambda).\n")
    }
    if(LT$type=="gamma")
    {
		if(is.null(LT$shape))
                {
			 LT$shape=1.1
			 if(verbose)
			 {
			 	cat(paste("  shape parameter in LP ",j," was missing and was set to ",LT$shape,"\n",sep=""))
			 }
		}
		
		if(is.null(LT$rate))
                {
			 LT$rate=(LT$shape-1)/LT$lambda2
			 if(verbose)
			 {
			 	cat(paste("  rate parameter in LP ",j," was missing and was set to ",LT$rate,"\n",sep=""))
			 }
		}	
    }
    
    if(LT$type=="beta")
    {
                if(is.null(LT$probIn))
  		{
    			LT$probIn=0.5
			if(verbose)
			{
    				cat(paste("  probIn in LP ",j," was missing and was set to ",LT$probIn,"\n",sep=""))
			}
  		}

  		if(is.null(LT$counts))
  		{
    			LT$counts=2
			if(verbose)
			{
    				cat(paste("  Counts in LP ",j," was missing and was set to ",LT$counts,"\n",sep=""))
			}
  		} 

                LT$shape1=LT$probIn*LT$counts;
                LT$shape2=(1-LT$probIn)*LT$counts;

		if(is.null(LT$max))
		{
		    LT$max=10*LT$lambda
		    if(verbose)
		    {
		    	cat(paste("  max parameter in LP ",j," was missing and was set to ",LT$max,"\n",sep=""))
		    }
		}
    }

    #Objects to storing information for MCMC iterations

    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    
    tmp=((var(y,na.rm=TRUE)*R2/nLT)/(LT$MSx))
    LT$tau2=rep(tmp,LT$p)
    LT$post_tau2=0  
    LT$post_lambda=0
    
    fname=paste(saveAt,LT$Name,"_lambda.dat",sep="");
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }
    
    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")

    LT$X=as.vector(LT$X)

    return(LT)
}


#Reproducing kernel Hilbert spaces
#This function simply sets hyperparameters and prepares inputs 
#for Reproducing Kernel Hilbert Spaces. The algorithm used here is 
#Fully described in de los Campos et al (2010).

setLT.RKHS=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles,verbose)  
{

    #Checking inputs
    if(is.null(LT$V))
    {
        if(is.null(LT$K)) stop(paste(" Kernel for linear term ",j, " was not provided, specify it with list(K=?,model='RKHS'), where ? is the kernel matrix\n",sep=""))

	LT$K = as.matrix(LT$K)

        if(class(LT$K)!="matrix") stop(paste(" Kernel for linear term ",j, " should be a matrix, the kernel provided is of class ", class(LT$K),"\n",sep=" "))

        if(nrow(LT$K)!=ncol(LT$K)) stop(paste(" Kernel for linear term ",j, " is not a square matrix\n",sep=""))

	#This code was rewritten to speed up computations
        #T = diag(weights)   
        #LT$K = T %*% LT$K %*% T 
        
        #Weight kernels
	for(i in 1:nrow(LT$K))
        {
		#j can not be used as subindex because its value is overwritten
		for(m in i:ncol(LT$K))
                {    
				LT$K[i,m]=LT$K[i,m]*weights[i]*weights[m];
                                LT$K[m,i]=LT$K[i,m]
		}
	}
        tmp =eigen(LT$K)
        LT$V =tmp$vectors
        LT$d =tmp$values
	rm(tmp)
    }else{
	if(any(weights!=1))
        { 
		cat(paste(" Warning, in LT ",j," Eigen decomposition was provided, and the model involves weights. Note: You should have weighted the kernel before computing eigen(K).\n",sep="")) 
        }
    }
    
    #Defaul value for tolD
    #Only those eigenvectors whose eigenvalues> tolD are kept.
    if (is.null(LT$tolD)) 
    {
       LT$tolD = 1e-10
       if(verbose)
       {
       		cat(paste("  Default value of minimum eigenvalue in LP ",j," was set to ",LT$tolD,"\n",sep=""))
       }
    }
    
    #Removing elements whose eigenvalues < tolD
    tmp= LT$d > LT$tolD
    LT$levelsU = sum(tmp)
    LT$d = LT$d[tmp]
    LT$V = LT$V[, tmp]
    
    #Default degrees of freedom and scale parameter associated with the variance component for marker effect
    if (is.null(LT$df0)) 
    {
      LT$df0 = 5
      if(verbose)
      {
      	cat(paste("  default value of df0 in LP ",j," was missing and was set to ",LT$df0,"\n",sep=""))
      }
    }
   
    if(is.null(LT$R2))
    { 
           LT$R2=R2/nLT
    }

    if (is.null(LT$S0)) 
    {
          if(LT$df0<=0) stop("df0>0 in RKHS in order to set S0\n");

	  LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(mean(LT$d)))*(LT$df0+2)

	  if(verbose)
	  {
             cat(paste("  default value of S0 in LP ",j," was missing and was set to ",LT$S0,"\n",sep=""))
	  }
    }
    
    LT$u=rep(0,nrow(LT$V))
    
    LT$varU=LT$S0/(LT$df0+2)
       
    LT$uStar=rep(0, LT$levelsU)
    
    #Output files
    fname=paste(saveAt,LT$Name,"_varU.dat",sep="")
    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    #Objects for storing information for MCMC iterations
    LT$fileOut=file(description=fname,open="w")
    LT$post_varU=0
    LT$post_varU2=0
    LT$post_uStar = rep(0, LT$levelsU)
    LT$post_u = rep(0, nrow(LT$V))
    LT$post_u2 = rep(0,nrow(LT$V))
    
    #return object
    return(LT)
}

###Bayes B and C########################################################################################################################################                 
#Pseudo BayesB with random scale and random proportion of markers "in" the model
#See Variable selection for regression models, 
#Lynn Kuo and Bani Mallic, 1998. 
#Bayes C (Habier et al., 2011)

setLT.BayesBandC=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles, groups, nGroups,verbose, nIter, burnIn, thin)
{

  model=LT$model
 
  if(is.null(LT$X)) LT$X=set.X(LT)

  #Be sure that your X is a matrix
  LT$X=as.matrix(LT$X)  
  LT$p=ncol(LT$X)
  LT$colNames=colnames(LT$X)

  #Weight inputs if necessary
  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights

  if(!is.null(groups))
  {

        x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
        for(g in 1:nGroups)
        {
                x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
        }
        LT$x2=x2;
    }else{
        LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  }

  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  
  if(any(is.na(LT$X))){ stop(paste("LP ",j," has NAs in X",sep=""))}
  if(nrow(LT$X)!=n){stop(paste("   Number of rows of LP ",j,"  not equal to the number of phenotypes.",sep=""))}
  
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    if(verbose)
    {
      cat(paste("  R2 in LP ",j," was missing and was set to ",LT$R2,"\n",sep=""))
    }
  }

  #Default value for the degrees of freedom associated with the distribution assigned to the variance
  #of marker effects
  if(is.null(LT$df0))
  {
    LT$df0= 5
    if(verbose)
    {
    	cat(paste("  DF in LP ",j," was missing and was set to ",LT$df0,"\n",sep=""))
    }
  }


  #Default value for a marker being "in" the model
  if(is.null(LT$probIn))
  {
    LT$probIn=0.5
    if(verbose)
    {	
       cat(paste("  probIn in LP ",j," was missing and was set to ",LT$probIn,"\n",sep=""))
    }
  } 


   #Default value for prior counts
  if(is.null(LT$counts))
  {
    LT$counts=10
    if(verbose)
    {
       cat(paste("  Counts in LP ",j," was missing and was set to ",LT$counts,"\n",sep=""))
    }
  }

  LT$countsIn=LT$counts * LT$probIn
  LT$countsOut=LT$counts - LT$countsIn

  #Default value for the scale parameter associated with the distribution assigned to the variance of 
  #marker effects
  if(is.null(LT$S0))
  {
     if(LT$df0<=0) stop(paste("df0>0 in ",model," in order to set S0\n",sep=""));
     LT$S0=var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0+2)/LT$probIn
     if(verbose)
     {
     	cat(paste(" Scale parameter in LP ",j," was missing and was set to ",LT$S0,"\n",sep=""))
     }
  }
 
  LT$b=rep(0, LT$p)
  LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)

  if(model=="BayesB")
  {
        if(is.null(LT$shape0))
  	{
        	LT$shape0=1.1
  	}
  	if(is.null(LT$rate0)){
    		LT$rate0=(LT$shape0-1)/LT$S0
  	}
        LT$S=LT$S0
  	LT$varB = LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
        fname=paste(saveAt,LT$Name,"_parBayesB.dat",sep="")
  }else{
	LT$varB = LT$S0
	fname=paste(saveAt,LT$Name,"_parBayesC.dat",sep="")
        
  }

  LT$X=as.vector(LT$X)
  LT$x2=as.vector(LT$x2)

  if(rmExistingFiles)
  { 
      unlink(fname) 
  }

  LT$fileOut=file(description=fname,open="w")
  LT$NamefileOut=fname;
  
  if(model=="BayesB")
  {
	tmp=c('probIn','scale')
   	write(tmp, ncolumns = LT$p, file = LT$fileOut, append = TRUE)
  }

  #Objects for storing MCMC information 
  LT$post_varB=0
  LT$post_varB2=0
  LT$post_d=0
  LT$MCMCinD = rep(0,LT$p)
  LT$BFw = rep(0,LT$p)
  LT$BF = rep(0,LT$p)
  LT$Xij = big.matrix(init=0, ncol=LT$p, nrow=((nIter-burnIn)/thin), type="char")
  LT$Pi = 0
  LT$post_probIn=0
  LT$post_probIn2=0
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)

  if(model=="BayesB")
  {
     LT$post_S=0
     LT$post_S2=0
  }
  
  #return object
  return(LT) 
}

#Bayes A, Mewissen et al. (2001).
#Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps
#Genetics 157: 1819-1829, Modified so that the Scale parameter is estimated from data (a gamma prior is assigned)

setLT.BayesA=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles,verbose)
{

  #Ckecking inputs 
  if(is.null(LT$X)) LT$X=set.X(LT)
 
  LT$X=as.matrix(LT$X)
  LT$p=ncol(LT$X)
  LT$colNames=colnames(LT$X)

  #Weight inputs if necessary
  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
  LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  #Default degrees of freedom for the prior assigned to the variance of markers
  if(is.null(LT$df0))
  {
     LT$df0 = 5 
     if(verbose)
     { 
     	cat(paste("  DF in LP ",j," was missing and was set to ",LT$df0,".\n",sep=""))
     }
  }
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    if(verbose)
    {
    	cat(paste("  R2 in LP ",j," was missing and was set to ",LT$R2,"\n",sep=""))
    }
  }

  #Defuault scale parameter for the prior assigned to the variance of markers
  if(is.null(LT$S0))
  {
     if(LT$df0<=0) stop("df0>0 in BayesA in order to set S0\n")
     LT$S0 = var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0+2)
     if(verbose)
     {
     	cat(paste(" Scale parameter in LP ",j," was missing and was set to ",LT$S0,"\n",sep=""))
     }
  }

  # Improvement: Treat Scale as random, assign a gamma density 
  if(is.null(LT$shape0))
  {
     LT$shape0=1.1
  }

  if(is.null(LT$rate0))
  {    
     LT$rate0=(LT$shape0-1)/LT$S0
  }
  LT$S=LT$S0
  
  LT$b=rep(0,LT$p)
   
  LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
  
  # Add one file when S0 is treated as random.
  fname=paste(saveAt,LT$Name,"_ScaleBayesA.dat",sep="") 
  if(rmExistingFiles)
  { 
    unlink(fname) 
  }

  LT$fileOut=file(description=fname,open="w")
  LT$NamefileOut=fname;
    
  LT$X=as.vector(LT$X)
  
  #Objects for storing information generated during MCMC iterations
  LT$post_varB=0
  LT$post_varB2=0
  
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)
  LT$post_S=0
  LT$post_S2=0
  
  #return object
  return(LT)
}

##################################################################################################
#Just the welcome function that will appear every time that your run the program
welcome=function()
{
  cat("\n");
  cat("#--------------------------------------------------------------------#\n");
  cat("#        _\\\\|//_                                                     #\n");
  cat("#       (` o-o ')      BGLR v1.0.4 build 98                          #\n");
  cat("#------ooO-(_)-Ooo---------------------------------------------------#\n");
  cat("#                      Bayesian Generalized Linear Regression        #\n");
  cat("#                      Gustavo de los Campos, gdeloscampos@gmail.com #\n");
  cat("#    .oooO     Oooo.   Paulino Perez, perpdgo@gmail.com              #\n");
  cat("#    (   )     (   )   April, 2015                                   #\n");
  cat("#_____\\ (_______) /_________________________________________________ #\n");
  cat("#      \\_)     (_/                                                   #\n");
  cat("#                                                                    #\n");
  cat("#------------------------------------------------------------------- #\n");
  cat("\n");
}
##################################################################################################

##################################################################################################
#The density of a scaled inverted chi-squered distribution
#df: degrees of freedom, S: Scale parameter
dScaledInvChisq=function (x, df, S)
{
    tmp = dchisq(S/x, df = df)/(x^2)
    return(tmp)
}


##################################################################################################
#The density function for lambda
#Density function for Regularization parameter in Bayesian LASSO
#Rate: rate parameter, shape: the value for the shape parameter
dLambda=function (rate, shape, lambda) 
{
    tmp = dgamma(x = I(lambda^2), rate = rate, shape = shape) * 2 * lambda
    return(tmp)
}

##################################################################################################
#Metropolis sampler for lambda in the Bayesian LASSO
metropLambda=function (tau2, lambda, shape1 = 1.2, shape2 = 1.2, max = 200, ncp = 0)
{
    lambda2 = lambda^2
    l2_new = rgamma(rate = sum(tau2)/2, shape = length(tau2),
        n = 1)
    l_new = sqrt(l2_new)
    logP_old = sum(dexp(x = tau2, log = TRUE, rate = (lambda2/2))) +
        dbeta(x = lambda/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/lambda2),
            log = TRUE)
    logP_new = sum(dexp(x = tau2, log = TRUE, rate = (l2_new/2))) +
        dbeta(x = l_new/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/l2_new),
            log = TRUE)
    accept = (logP_new - logP_old) > log(runif(1))
    if (accept) {
        lambda = l_new
    }
    return(lambda)
}

##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.1.3"))
    stop("This package requires R 3.1.3 or later")
  assign(".BGLR.home", file.path(library, pkg),
         pos=match("package:BGLR", search()))
  BGLR.version = "1.0.4 (2015-04-07), build 98"
  assign(".BGLR.version", BGLR.version, pos=match("package:BGLR", search()))
  if(interactive())
  {
    packageStartupMessage(paste("# Package Bayesian Generalized Regression (BGLR), ", BGLR.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("# Gustavo de los Campos & Paulino Perez",appendLF=TRUE)
    packageStartupMessage("# Support provided by the U.S., National Institutes of Health (NIH)", appendLF=TRUE)
    packageStartupMessage("# (Grant: R01GM101219, NIGMS)", appendLF=TRUE)
    packageStartupMessage("# and by the International Maize and Wheat Improvement Center (CIMMyT).",appendLF=TRUE)                        
    packageStartupMessage("# Type 'help(BGLR)' for summary information",appendLF=TRUE)
  }
  invisible()
}


##################################################################################################
#rtrun draws from a truncated univariate normal distribution using the inverse CDF algorithm
#Arguments:
#mu: mean
#sigma: standard deviation
#a: lower bound
#b: upper bound
#NOTES: 1) This routine was taken from bayesm package, December 18, 2012
#       2) The inputs are not checked, 
#It is assumed that are ok.
rtrun=function (mu, sigma, a, b) 
{
    FA = pnorm(((a - mu)/sigma))
    FB = pnorm(((b - mu)/sigma))
    return(mu + sigma * qnorm(runif(length(mu)) * (FB - FA) + FA))
}

#Extract the values of z such that y[i]=j
#z,y vectors, j integer
extract=function(z,y,j) subset(as.data.frame(z,y),subset=(y==j))

#This routine was adapted from rinvGauss function from S-Plus
# Random variates from inverse Gaussian distribution
# Reference:
#      Chhikara and Folks, The Inverse Gaussian Distribution,
#      Marcel Dekker, 1989, page 53.
# GKS  15 Jan 98
#n: Number of samples
#nu: nu parameter
#lambda: lambda parameter

rinvGauss=function(n, nu, lambda)
{
	if(any(nu<=0)) stop("nu must be positive")
	if(any(lambda<=0)) stop("lambda must be positive")
	if(length(n)>1) n = length(n)
	if(length(nu)>1 && length(nu)!=n) nu = rep(nu,length=n)
	if(length(lambda)>1 && length(lambda)!=n) lambda = rep(lambda,length=n)
        tmp = rnorm(n)
	y2 = tmp*tmp
	u = runif(n)
	r1 = nu/(2*lambda) * (2*lambda + nu*y2 - sqrt(4*lambda*nu*y2 + nu*nu*y2*y2))
	r2 = nu*nu/r1
	ifelse(u < nu/(nu+r1), r1, r2)
}


#log-likelihood for ordinal data
#y: response vector
#predicted response vector, yHat=X%*%beta
#threshold
loglik_ordinal=function(y,yHat,threshold)
{
	sum=0
        n=length(y)
	for(i in 1:n)
        {
           sum=sum + log(pnorm(threshold[y[i] + 1]-yHat[i])-pnorm(threshold[y[i]]-yHat[i]))
        }
        return(sum)
}

##################################################################################################

#Arguments:
#y: data vector, NAs allowed
#response_type: can be "gaussian", "ordinal",
#ETA: The linear predictor
#nIter: Number of MCMC iterations
#burnIn: burnIn
#thin: thin
#saveAt: string, where to save the information
#Se: Scale parameter for the prior for varE
#dfe: Degrees of freedom for the prior for varE
#weights: 
#R2
#Note: The function was designed to work with gaussian responses, some changes were made to deal binary and ordinal responses


#To add new method: 
#(a) create setLT, 
#(b) add it to the switch statement,
#(c) add code to update parameters in the Gibbs sampler, 
#(d) add code to save samples
#(e) add code to compute posterior means
#(f) Test:
#(f1) Test simple example without hyper-paramaeters, evaluate how  
#        default values were set
	#(f2)  Check posterior means and files
	#(f3)  Test simple example with some hyper-parameters give and 
	#         some set by default
#(f4) Run an example with a few missing values, compare to BLR 
#       example, check: (i) residual variance, (ii) plot of effects, (iii) plot 
#        of predictions in trn, (iv) plot of prediction in tst.


BGLR=function (y, response_type = "gaussian", a = NULL, b = NULL, 
    ETA = NULL, nIter = 1500, burnIn = 500, thin = 1, saveAt = "", 
    S0 = NULL, df0 = 5, R2 = 0.5, weights = NULL, 
    verbose = TRUE, rmExistingFiles = TRUE, groups=NULL, BigMatrixName=NULL, map=map, CustomLimit=4, Windowed=TRUE) 
{
   
    if(verbose)
    {
	welcome()
    }
	
	if("bigmemory" %in% rownames(installed.packages())==FALSE){
		cat("Script requires package 'bigmemory' ", "\n")
		n<-readline("Should the package be installed now? (If no then the function will abort.) (y/n): ")
		if(n=="y"){
			install.packages("bigmemory")
		} else {
			return()
		}
	}
	if("bigmemory" %in% rownames(search())==FALSE){
		library(bigmemory)
	}
	options(bigmemory.typecast.warning=FALSE)
    IDs=names(y)
    if (!(response_type %in% c("gaussian", "ordinal")))  stop(" Only gaussian and ordinal responses are allowed\n")

    if (saveAt == "") {
        saveAt = paste(getwd(), "/", sep = "")
    }

    y=as.vector(y)
    y0=y
    a = as.vector(a)
    b = as.vector(b)
    n = length(y)

    nGroups=1
    if(!is.null(groups))
    {
		groups<-as.character(groups)  #Groups as character and then as factor to avoid dummy levels
		groups<-as.factor(groups)
		#Number of records by group
		countGroups=table(groups)
		nGroups=length(countGroups)
		groupLabels=names(countGroups)
                groups=as.integer(groups)
                ggg=as.integer(groups-1);  #In C we begin to count in 0
		if(sum(countGroups)!=n) stop("length of groups and y differs, NA's not allowed in groups\n");	
    }

    if(response_type=="ordinal")
    {

    	y=factor(y,ordered=TRUE)
        lev=levels(y)
        nclass=length(lev)
        if(nclass==n) stop("The number of classes in y must be smaller than the number of observations\n");

        y=as.integer(y)
        z=y  
    }
    
    if (is.null(weights)) 
    {
        weights = rep(1, n)
    }

    if(!is.null(groups))
    {
      sumW2=tapply(weights^2,groups,"sum")
    }else{
      sumW2 = sum(weights^2)
    }

    nSums = 0
	Pival=0.5 #Only used in BF calculation. Updated by probIn.
    whichNa = which(is.na(y))
    nNa = length(whichNa)
	IOvector = rep() 
	
    Censored = FALSE
     
    if (response_type == "gaussian") 
    {
        if ((!is.null(a)) | (!is.null(b))) 
        {
            Censored = TRUE
            if ((length(a) != n) | (length(b) != n)) stop(" y, a and b must have the same dimension\n")
            if (any(weights != 1)) stop(" Weights are only implemented for Gausian uncensored responses\n")
        }
        mu = weighted.mean(x = y, w = weights, na.rm = TRUE)
    }
    post_mu = 0
    post_mu2 = 0

    fname = paste(saveAt, "mu.dat", sep = "")
    if (rmExistingFiles) 
    {
        unlink(fname)
    }
    else {
        cat(" Note: samples will be appended to existing files. \n")
    }

    fileOutMu = file(description = fname, open = "w")

    if (response_type == "ordinal") {
        cat(" Prior for residual is not necessary, if you provided it, it will be ignored\n")
        if (any(weights != 1)) stop(" Weights are not supported \n")
       
        countsZ=table(z)

        if (nclass <= 1) stop(paste(" Data vector y has only ", nclass, " differente values, it should have at least 2 different values\n"))
        threshold=qnorm(p=c(0,cumsum(as.vector(countsZ)/n)))
          
        y = rtrun(mu =0, sigma = 1, a = threshold[z], b = threshold[ (z + 1)])
        mu=0
        #posterior for thresholds
        post_threshold = 0
        post_threshold2 = 0
        
	post_prob=matrix(nrow=n,ncol=nclass,0)
        post_prob2=post_prob
    }

    post_logLik = 0

    # yStar & yHat
    yStar = y * weights
    yHat = mu * weights
    
    if (nNa > 0) {
        yStar[whichNa] = yHat[whichNa]
    }

    post_yHat = rep(0, n)
    post_yHat2 = rep(0, n)

    # residual and residual variance
    e = (yStar - yHat)

    varE = var(e, na.rm = TRUE) * (1 - R2)

    if (is.null(S0)) {
        S0 = varE * (df0 + 2)
    }

    if(!is.null(groups))
    {
        varE=rep(varE/nGroups,nGroups)
        names(varE)=groupLabels
    }

    sdE = sqrt(varE)


    post_varE = 0
    post_varE2 = 0

    #File for storing sample for varE 

    fname = paste(saveAt, "varE.dat", sep = "")

    if (rmExistingFiles) {
        unlink(fname)
    }

    fileOutVarE = file(description = fname, open = "w")

    nLT = ifelse(is.null(ETA), 0, length(ETA))
    

    #Setting the linear terms
    if (nLT > 0) {
	
	if(is.null(names(ETA)))
    	{ 
             names(ETA)<-rep("",nLT)
    	}

        for (i in 1:nLT) {  

	    if(names(ETA)[i]=="")
	    {
	       	ETA[[i]]$Name=paste("ETA_",i,sep="")
	    }else{
               ETA[[i]]$Name=paste("ETA_",names(ETA)[i],sep="")
	    }

            if (!(ETA[[i]]$model %in% c("FIXED", "BRR", "BL", "BayesA", "BayesB","BayesC", "RKHS","BRR_windows"))) 
            {
                stop(paste(" Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented (note: evaluation is case sensitive).", sep = ""))
                
            }

            if(!is.null(groups))
            {
		if(!(ETA[[i]]$model %in%  c("BRR","FIXED","BayesB","BayesC"))) stop(paste(" Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented for groups\n", sep = ""))
            }


            ETA[[i]] = switch(ETA[[i]]$model, 
			      FIXED = setLT.Fixed(LT = ETA[[i]],  n = n, j = i, weights = weights, y = y, nLT = nLT, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups), 
                              BRR = setLT.BRR(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose), 
                              BL = setLT.BL(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose), 
                              RKHS = setLT.RKHS(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose), 
                              BayesC = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose, nIter=nIter, burnIn=burnIn, thin=thin),
                              BayesA = setLT.BayesA(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose),
                              BayesB = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose, nIter=nIter, burnIn=burnIn, thin=thin),
                              BRR_windows = setLT.BRR_windows(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose)
                              )
        }
    }

    # Gibbs sampler

    time = proc.time()[3]

    for (i in 1:nIter) {
        # intercept
	if(!is.null(groups))
	{
		e = e + weights * mu 
                varEexpanded=varE[groups]
		#rhs = sum(tapply(e*weights,groups,"sum")/varE)
                rhs = as.numeric(crossprod(e/varEexpanded,weights));
                C = sum(sumW2/varE)
                sol = rhs/C
                mu = rnorm(n = 1, sd = sqrt(1/C)) + sol;
	}else{
        	e = e + weights * mu
        	rhs = sum(weights * e)/varE
        	C = sumW2/varE
        	sol = rhs/C
        	mu = rnorm(n = 1, sd = sqrt(1/C)) + sol
	}
        if (response_type == "ordinal") {
            mu=0 
        }

        e = e - weights * mu #residuaali-keskiarvo?
        
        #deltaSS and deltadf for updating varE
        deltaSS = 0
        deltadf = 0

        if (nLT > 0) {
            for (j in 1:nLT) {
                ## Fixed effects ####################################################################
                if (ETA[[j]]$model == "FIXED") {
                  #cat("varB=",ETA[[j]]$varB,"\n");
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)
                  if(!is.null(groups)){
                        ans = .Call("sample_beta_groups", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                             e, varBj, varE, 1e-9,ggg,nGroups)
		  }else{
                  	ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, 1e-9)
		  }
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                }#End of fixed effects

                ## Ridge Regression ##################################################################
                if (ETA[[j]]$model == "BRR") {
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)

                  if(!is.null(groups))
		  {
                        ans = .Call("sample_beta_groups",n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                    e, varBj, varE, 1e-9,ggg,nGroups)
	          }else{
                  	ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, 1e-9)
		  }
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  DF = ETA[[j]]$df0 + ETA[[j]]$p
                  SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
                  ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
                }# END BRR
                
                if(ETA[[j]]$model=="BRR_windows"){
                   ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                             e, ETA[[j]]$varB, varE, 1e-9)
		   ETA[[j]]$b = ans[[1]]
                   e = ans[[2]]

                   tmp=numeric()
                   for(nw in 1:ETA[[j]]$nwindows)
                   {
                        index=ETA[[j]]$windows_list[[nw]]
			DF = ETA[[j]]$df0 + length(index)
                  	SS = sum((ETA[[j]]$b[index])^2) + ETA[[j]]$S0
                        tmp=c(tmp,rep(SS/rchisq(df = DF, n = 1),length(index)))
                   }
                   ETA[[j]]$varB=tmp
                }


                ## Bayesian LASSO ####################################################################
                if (ETA[[j]]$model == "BL") {
                  
                   varBj = ETA[[j]]$tau2 * varE
                   ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                e, varBj, varE, ETA[[j]]$minAbsBeta)

                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  nu = sqrt(varE) * ETA[[j]]$lambda/abs(ETA[[j]]$b)
                  tmp = NULL
                  try(tmp <- rinvGauss(n = ETA[[j]]$p, nu = nu, lambda = ETA[[j]]$lambda2))
                  if (!is.null(tmp) && !any(tmp<0)) {
                    if (!any(is.na(sqrt(tmp)))) {
                      ETA[[j]]$tau2 = 1/tmp
                    }
                    else {
                      warning(paste("tau2 was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }
                  else {
                    warning(paste("tau2 was not updated  in iteration",i,"due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                  }

                  #Update lambda 
                  if (ETA[[j]]$type == "gamma") {
                    rate = sum(ETA[[j]]$tau2)/2 + ETA[[j]]$rate
                    shape = ETA[[j]]$p + ETA[[j]]$shape
                    ETA[[j]]$lambda2 = rgamma(rate = rate, shape = shape, n = 1)
                    if (!is.na(ETA[[j]]$lambda2)) {
                      ETA[[j]]$lambda = sqrt(ETA[[j]]$lambda2)
                    }
                    else {
                      warning(paste("lambda was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }

                  if (ETA[[j]]$type == "beta") {
                    ETA[[j]]$lambda = metropLambda(tau2 = ETA[[j]]$tau2, 
                                                   lambda = ETA[[j]]$lambda, shape1 = ETA[[j]]$shape1, shape2 = ETA[[j]]$shape2, 
                                                   max = ETA[[j]]$max)
                    ETA[[j]]$lambda2 = ETA[[j]]$lambda^2
                  }

                  deltaSS = deltaSS + sum((ETA[[j]]$b/sqrt(ETA[[j]]$tau2))^2)
                  deltadf = deltadf + ETA[[j]]$p
                }#END BL

                ## RKHS ####################################################################
                if (ETA[[j]]$model == "RKHS") {
                  #error
                  e = e + ETA[[j]]$u
                  rhs = crossprod(ETA[[j]]$V, e)/varE
                  varU = ETA[[j]]$varU * ETA[[j]]$d
                  C = as.numeric(1/varU + 1/varE)
                  SD = 1/sqrt(C)
                  sol = rhs/C
                  tmp = rnorm(n = ETA[[j]]$levelsU, mean = sol, sd = SD)
                  ETA[[j]]$uStar = tmp
                  ETA[[j]]$u = as.vector(ETA[[j]]$V %*% tmp)
		  
                  #update error
                  e = e - ETA[[j]]$u
                   
                  #update the variance
                  tmp = ETA[[j]]$uStar/sqrt(ETA[[j]]$d)
                  SS = as.numeric(crossprod(tmp)) + ETA[[j]]$S0
                  DF = ETA[[j]]$levelsU + ETA[[j]]$df0
                  ETA[[j]]$varU = SS/rchisq(n = 1, df = DF)
                }#END RKHS

                ## BayesA ##############################################################################
                if (ETA[[j]]$model == "BayesA") {
                  varBj = ETA[[j]]$varB
                  ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, 1e-9)
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                   
                  #Update variances
                  SS = ETA[[j]]$S + ETA[[j]]$b^2
                  DF = ETA[[j]]$df0 + 1
                  ETA[[j]]$varB = SS/rchisq(n = ETA[[j]]$p, df = DF)

                  tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                  tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                  ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)

                }#End BayesA

		#BayesB and BayesC
		if(ETA[[j]]$model %in% c("BayesB","BayesC"))
		{
			Pival = ETA[[j]]$probIn
			#Update marker effects
                      	mrkIn=ETA[[j]]$d==1
                      	pIn=sum(mrkIn)
		        
                        if(ETA[[j]]$model=="BayesB")
                        {
                          if(!is.null(groups))
                          {
                             ans=.Call("sample_beta_BB_BCp_groups",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, ETA[[j]]$varB, varE, 1e-9, ETA[[j]]$probIn,ggg,nGroups);
                          }else{
                             ans=.Call("sample_beta_BB_BCp",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, ETA[[j]]$varB, varE, 1e-9, ETA[[j]]$probIn);
                          }
                        }else{
                          if(!is.null(groups))
                          {
                             ans=.Call("sample_beta_BB_BCp_groups",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, rep(ETA[[j]]$varB,ETA[[j]]$p), varE, 1e-9, ETA[[j]]$probIn,ggg,nGroups);
                          }else{   
                             ans=.Call("sample_beta_BB_BCp",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, rep(ETA[[j]]$varB,ETA[[j]]$p), varE, 1e-9, ETA[[j]]$probIn);
                          }
                        }

                        ETA[[j]]$d=ans[[1]]
                        e=ans[[2]]
                        ETA[[j]]$b=ans[[3]] 

			#Update the variance component associated with the markers
			if(ETA[[j]]$model=="BayesB")
			{
				SS = ETA[[j]]$b^2 + ETA[[j]]$S
                      		DF = ETA[[j]]$df0+1
                      		ETA[[j]]$varB = SS/rchisq(df=DF, n = ETA[[j]]$p)

                                # Update scale
                                tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                                tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                                ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)

			}else{
				SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
				DF = ETA[[j]]$df0 + ETA[[j]]$p
                  		ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
			}
                      	mrkIn = sum(ETA[[j]]$d)
                      	ETA[[j]]$probIn = rbeta(shape1 = (mrkIn + ETA[[j]]$countsIn + 1),
                                                shape2 = (ETA[[j]]$p - mrkIn + ETA[[j]]$countsOut + 1), n = 1)
						
						
		}
            }#Loop for
        }#nLT
        
        # yHat
        yHat = yStar - e
        
       #4#
         # residual variance # missing values
        if (response_type == "gaussian") {
	    
            if(!is.null(groups))
	    {	
        	for(g in 1:nGroups)
        	{
			SS=sum(e[groups==g]^2)+ S0 + deltaSS
                        DF=countGroups[g]+df0+deltadf
			varE[g]=SS/rchisq(n=1,df=DF)
		}
	     }else{
            		SS = sum(e * e) + S0 + deltaSS
            		DF = n + df0 + deltadf
            		varE = SS/rchisq(n = 1, df = DF)
	    }
            sdE = sqrt(varE)
        
          if (nNa > 0) {
            if (Censored) {
                if(!is.null(groups))
                {
                   #FIXME: Double check this, I was testing it and is ok
                   sdEexpanded=sdE[groups]
                   yStar[whichNa] = rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdEexpanded)

                }else{
                  yStar[whichNa] = rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdE)
                }
            }
            else{
                 if(!is.null(groups))
                 {
                    #FIXME: Double check this, I was testing it and is ok
                    sdEexpanded=sdE[groups]
                    yStar[whichNa] = yHat[whichNa] + rnorm(n = nNa, sd = sdEexpanded)
                 }else{
                    yStar[whichNa] = yHat[whichNa] + rnorm(n = nNa, sd = sdE)
                 }
            }
            e[whichNa] = yStar[whichNa] - yHat[whichNa]
          }
        }else{  #ordinal
            varE = 1
            sdE = 1
            
            #Update yStar, this is the latent variable
            if(nNa==0){
               yStar=rtrun(mu = yHat, sigma = 1, a = threshold[z], b = threshold[(z + 1)])
            }else{
               yStar[-whichNa]=rtrun(mu = yHat[-whichNa], sigma = 1, a = threshold[z[-whichNa]], b = threshold[(z[-whichNa] + 1)])
               yStar[whichNa]=yHat[whichNa] + rnorm(n = nNa, sd = sdE)           
            }

            #Update thresholds           
            if(nNa==0){ 
              for (m in 2:nclass) {
            
                lo = max(max(extract(yStar, z, m - 1)), threshold[m - 1])
                hi = min(min(extract(yStar, z, m)), threshold[m + 1])
                threshold[m] = runif(1, lo, hi)
              }
            }else{

              for (m in 2:nclass) {
                tmpY=yStar[-whichNa]
                tmpZ=z[-whichNa]
                lo = max(max(extract(tmpY, tmpZ, m - 1)), threshold[m - 1])
                hi = min(min(extract(tmpY, tmpZ, m)), threshold[m + 1])
                threshold[m] = runif(1, lo, hi)
              }
            }
            
            #Update error
            e = yStar - yHat
        }

        # Saving samples and computing running means
        if ((i%%thin == 0)) {
            if (nLT > 0) {
                for (j in 1:nLT) {

                  if (ETA[[j]]$model == "FIXED") {
                    write(ETA[[j]]$b,ncolumns=ETA[[j]]$p, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BRR") {
                    write(ETA[[j]]$varB, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BL") {
                    write(ETA[[j]]$lambda, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "RKHS") {
                    write(ETA[[j]]$varU, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesC") {
                    tmp = c(ETA[[j]]$probIn, ETA[[j]]$varB)
                    write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesA") {
                    tmp=ETA[[j]]$S
                    write(tmp, ncolumns = 1, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                  
                  if(ETA[[j]]$model=="BayesB")
                  {
                        tmp=c(ETA[[j]]$probIn,ETA[[j]]$S)
                        write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                }
            }

            #Output files
            write(x = mu, file = fileOutMu, append = TRUE)
            write(x = varE, ncolumns=nGroups,file = fileOutVarE, append = TRUE)
            if (i > burnIn) {
                nSums = nSums + 1
                k = (nSums - 1)/(nSums)
                if (nLT > 0) {
                  for (j in 1:nLT) {

                    if (ETA[[j]]$model == "FIXED") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                    }

                    if (ETA[[j]]$model == "BRR") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                    }

                    if (ETA[[j]]$model == "BL") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_tau2 = ETA[[j]]$post_tau2 * k + (ETA[[j]]$tau2)/nSums
                      ETA[[j]]$post_lambda = ETA[[j]]$post_lambda * k + (ETA[[j]]$lambda)/nSums
                    }

                    if (ETA[[j]]$model == "RKHS") {
                      ETA[[j]]$post_varU = ETA[[j]]$post_varU * k + ETA[[j]]$varU/nSums
                      ETA[[j]]$post_varU2 = ETA[[j]]$post_varU2 * k + (ETA[[j]]$varU^2)/nSums
                      ETA[[j]]$post_uStar = ETA[[j]]$post_uStar * k + ETA[[j]]$uStar/nSums
                      ETA[[j]]$post_u = ETA[[j]]$post_u * k + ETA[[j]]$u/nSums
                      ETA[[j]]$post_u2 = ETA[[j]]$post_u2 * k + (ETA[[j]]$u^2)/nSums
                    }

                    if (ETA[[j]]$model == "BayesC") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$Pi <- append(ETA[[j]]$Pi,Pival)
					  ETA[[j]]$MCMCinD = ETA[[j]]$MCMCinD + ETA[[j]]$d
					  ETA[[j]]$BF = ETA[[j]]$BF + ((ETA[[j]]$d)/(Pival/(1-Pival)))
					  ETA[[j]]$Xij[((i-burnIn)/thin), ETA[[j]]$d==1] <- 1
					  if(sum(ETA[[j]]$Xij[((i-burnIn)/thin),])==0){
							cat("Writing to big.matrix object failed at iteration", i,"\n", "returning", "\n")
							return()
						}
					  ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                      ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                      ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums
                    }

                    if (ETA[[j]]$model == "BayesA") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_S = ETA[[j]]$post_S * k + (ETA[[j]]$S)/nSums
					  ETA[[j]]$post_S2 = ETA[[j]]$post_S2 * k + (ETA[[j]]$S^2)/nSums
                    }

                    if(ETA[[j]]$model=="BayesB")
                    {
                        ETA[[j]]$post_b=ETA[[j]]$post_b*k+ETA[[j]]$b/nSums
                        ETA[[j]]$post_b2=ETA[[j]]$post_b2*k+(ETA[[j]]$b^2)/nSums
                        ETA[[j]]$post_varB=ETA[[j]]$post_varB*k+(ETA[[j]]$varB)/nSums
                        ETA[[j]]$post_varB2=ETA[[j]]$post_varB2*k+(ETA[[j]]$varB^2)/nSums
                        ETA[[j]]$Pi <- append(ETA[[j]]$Pi,Pival)
						ETA[[j]]$MCMCinD = ETA[[j]]$MCMCinD + ETA[[j]]$d
						ETA[[j]]$BF = ETA[[j]]$BF + ((ETA[[j]]$d)/(Pival/(1-Pival)))
						ETA[[j]]$Xij[((i-burnIn)/thin), ETA[[j]]$d==1] <- 1
						if(sum(ETA[[j]]$Xij[((i-burnIn)/thin),])==0){
							cat("Writing to big.matrix object failed at iteration", i,"\n", "returning", "\n")
							return()
						}
						ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                        ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                        ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums
                        ETA[[j]]$post_S = ETA[[j]]$post_S * k + (ETA[[j]]$S)/nSums
						ETA[[j]]$post_S2 = ETA[[j]]$post_S2 * k + (ETA[[j]]$S^2)/nSums
                    }
                  }
                }

                post_mu = post_mu * k + mu/nSums
                post_mu2 = post_mu2 * k + (mu^2)/nSums

                post_yHat = post_yHat * k + yHat/nSums
                post_yHat2 = post_yHat2 * k + (yHat^2)/nSums

                post_varE = post_varE * k + varE/nSums
                post_varE2 = post_varE2 * k + (varE^2)/nSums

                if (response_type == "ordinal") {
                  post_threshold = post_threshold * k + threshold/nSums
                  post_threshold2 = post_threshold2 * k + (threshold^2)/nSums

                  TMP=matrix(nrow=n,ncol=nclass,0)

                  TMP[,1]=pnorm(threshold[2]-yHat)

                  if(nclass>2){
                     for(m in 2:(nclass-1)){
                       TMP[,m]=pnorm(threshold[(m+1)]-yHat)-rowSums(as.matrix(TMP[,1:(m-1)]))
                     }
                  }
                  TMP[,nclass]=1-rowSums(TMP)

                  post_prob=post_prob*k+TMP/nSums
                  post_prob2=post_prob2*k+(TMP^2)/nSums
                 
                  if(nNa==0){
                    logLik=loglik_ordinal(z,yHat,threshold)
		          }else{
		            logLik=loglik_ordinal(z[-whichNa],yHat[-whichNa],threshold)
		          }
                }

                if(response_type == "gaussian") {
                  
                  tmpE = e/weights
                  if(!is.null(groups))
		  {
                    tmpSD=rep(NA,n)
                    for(g in 1:nGroups)
                    {
                       index=(groups==g)
                       tmpSD[index]=sqrt(varE[g])/weights[index]
                    }
                  }else{
                    tmpSD = sqrt(varE)/weights
		  }

                  if (nNa > 0) {
                    tmpE = tmpE[-whichNa]
                    tmpSD = tmpSD[-whichNa]
                  }
                  
	          logLik = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

                }#end gaussian

                post_logLik = post_logLik * k + logLik/nSums
            }
        }#end of saving samples and computing running means

        if (verbose) {
            cat("---------------------------------------\n")
            tmp = proc.time()[3]
            cat(c(paste(c("  Iter=", "Time/Iter="), round(c(i, c(tmp - time)), 3), sep = "")), "\n")
            cat("  VarE=",round(varE,3),"\n")
            time = tmp
        }
    }#end of Gibbs sampler

    #Closing files
    close(fileOutVarE)
    close(fileOutMu)

    if (nLT > 0) {
        for (i in 1:nLT) {
            if (!is.null(ETA[[i]]$fileOut)) {
                close(ETA[[i]]$fileOut)
            }
            ETA[[i]]$fileOut = NULL
        }
    }
    
    #return goodies
	
		
    out = list(y = y0, whichNa = whichNa, saveAt = saveAt, nIter = nIter, 
               burnIn = burnIn, thin = thin, 
               weights = weights, verbose = verbose, 
               response_type = response_type, df0 = df0, S0 = S0)

    out$yHat = post_yHat

    names(out$yHat)=IDs
    names(out$y)=IDs

    out$SD.yHat = sqrt(post_yHat2 - (post_yHat^2))
    out$mu = post_mu
    out$SD.mu = sqrt(post_mu2 - post_mu^2)
    out$varE = post_varE
    out$SD.varE = sqrt(post_varE2 - post_varE^2)
    
    #goodness of fit 
    out$fit = list()
    
    if(response_type=="gaussian")
    {
    	tmpE = (yStar - post_yHat)/weights

        if(!is.null(groups))
        {
                    tmpSD=rep(NA,n)
                    for(g in 1:nGroups)
                    {
                       index=(groups==g)
                       tmpSD[index]=sqrt(varE[g])/weights[index]
                    }
         }else{
    		tmpSD = sqrt(post_varE)/weights
	 }
    
    	if (nNa > 0) {
        	tmpE = tmpE[-whichNa]
        	tmpSD = tmpSD[-whichNa]
    	}
    	out$fit$logLikAtPostMean = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

	if (Censored) {
            cdfA = pnorm(q = a[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            cdfB = pnorm(q = b[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            out$fit$logLikAtPostMean = out$fit$logLikAtPostMean + sum(log(cdfB - cdfA))
        }
    }

    if(response_type=="ordinal")
    {
         out$fit$logLikAtPostMean = loglik_ordinal(y,post_yHat,post_threshold)
         out$probs=post_prob
         out$SD.probs=sqrt(post_prob2-post_prob^2)
         colnames(out$probs)=lev
         colnames(out$SD.probs)=lev
         out$threshold = post_threshold[-c(1, nclass + 1)]
         out$SD.threshold = sqrt(post_threshold2 - post_threshold^2)[-c(1, nclass + 1)] 

         out$levels=lev
         out$nlevels=nclass
    }

    out$fit$postMeanLogLik = post_logLik
    out$fit$pD = -2 * (post_logLik - out$fit$logLikAtPostMean)
    out$fit$DIC = out$fit$pD - 2 * post_logLik

    # Renaming/removing objects in ETA and appending names
    if (nLT > 0) {
        for (i in 1:nLT) {

            if (ETA[[i]]$model != "RKHS") {
                ETA[[i]]$b = ETA[[i]]$post_b
                ETA[[i]]$SD.b = sqrt(ETA[[i]]$post_b2 - ETA[[i]]$post_b^2)
                names(ETA[[i]]$b)=ETA[[i]]$colNames
                names(ETA[[i]]$SD.b)=ETA[[i]]$colNames
                tmp = which(names(ETA[[i]]) %in% c("post_b", "post_b2","X","x2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
            
            if(ETA[[i]]$model=="RKHS")
            {
               ETA[[i]]$SD.u=sqrt(ETA[[i]]$post_u2 - ETA[[i]]$post_u^2)
               ETA[[i]]$u=ETA[[i]]$post_u
               ETA[[i]]$uStar=ETA[[i]]$post_uStar
               ETA[[i]]$varU=ETA[[i]]$post_varU
               ETA[[i]]$SD.varU=sqrt(ETA[[i]]$post_varU2 - ETA[[i]]$post_varU^2)
               tmp=which(names(ETA[[i]])%in%c("post_varU","post_varU2","post_uStar","post_u","post_u2"))
               ETA[[i]]=ETA[[i]][-tmp]
            }

            if (ETA[[i]]$model %in% c("BRR", "BayesA", "BayesC","BayesB")) {
                ETA[[i]]$varB = ETA[[i]]$post_varB
                ETA[[i]]$SD.varB = sqrt(ETA[[i]]$post_varB2 - (ETA[[i]]$post_varB^2))
                tmp = which(names(ETA[[i]]) %in% c("post_varB", "post_varB2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }

            if(ETA[[i]]$model %in% c("BayesB","BayesC"))
            {
	        
				ETA[[i]]$BF = ETA[[i]]$BF * (1/(((nIter-burnIn)/thin)-ETA[[i]]$MCMCinD))
				if(is.null(BigMatrixName)==FALSE){
					write.big.matrix(x=ETA[[i]]$Xij, filename=BigMatrixName, row.names = FALSE, col.names = FALSE,sep = ",")
					cat("The indicator matrix is stored in ", BigMatrixName, ", defined in the 'BigMatrixName' in the modified BGLR function", "\n",
					"and can be read to R with the read.big.matrix function in the bigmemory R-package.", "\n", "(y-axis=iteration, x-axis=covariates, type=char)", "\n")
				}
			
				ETA[[i]]$Pi <- ETA[[i]]$Pi[-1] 
				ETA[[i]]$d=ETA[[i]]$post_d
                ETA[[i]]$probIn=ETA[[i]]$post_probIn
				ETA[[i]]$SD.probIn=sqrt(ETA[[i]]$post_probIn2 - (ETA[[i]]$post_probIn^2))
                tmp = which(names(ETA[[i]]) %in% c("post_d", "post_probIn","post_probIn2"))
                ETA[[i]] = ETA[[i]][-tmp]
				BF<-ETA[[i]]$BF
				b<-ETA[[i]]$b
				Pi<-ETA[[i]]$Pi
				MCMCinD<-ETA[[i]]$MCMCinD
				SignfLimKa=(sum(1/(Pi/(1-Pi))))/(nIter-burnIn)
				#SignfLimKa=1/((Pi[order(Pi)][length(Pi)*0.05])/(1-Pi[order(Pi)][length(Pi)*0.05]))
				SignfLim=1/((Pi[order(Pi)][1])/(1-Pi[order(Pi)][1]))
				ETA[[i]]$Xij<-getResults(indicatorMatrx=ETA[[i]]$Xij, map=map, BFlimit=SignfLim, BFlimit2=SignfLimKa, BF=BF, b=b, Pi=Pi, BFlimitCustom=CustomLimit)
            
			}
            
            if(ETA[[i]]$model %in% c("BayesA","BayesB"))
            {
                ETA[[i]]$S=ETA[[i]]$post_S
                ETA[[i]]$SD.S=sqrt( ETA[[i]]$post_S2 - (ETA[[i]]$post_S^2))
                tmp=which(names(ETA[[i]])%in%c("post_S","post_S2"))
                ETA[[i]]=ETA[[i]][-tmp]
            }

            if(ETA[[i]]$model=="BL")
            {
                ETA[[i]]$tau2=ETA[[i]]$post_tau2
                ETA[[i]]$lambda=ETA[[i]]$post_lambda
                tmp = which(names(ETA[[i]]) %in% c("post_tau2", "post_lambda","lambda2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
        }
        out$ETA = ETA
    }
	options(bigmemory.typecast.warning=TRUE)
    class(out) = "BGLR"
    return(out)
}

 #fm<-BGLR(y=y,ETA=list(list(X=X,model='FIXED')))

#This function will be a wrapper for BGLR
#the idea is to maintain the compatibility with the function BLR in 
#the package BLR that was released in 2010, updated in 2011 and 2012

#NOTE: thin2 parameter is missing in BGLR, so it will be removed

BLR=function (y, XF = NULL, XR = NULL, XL = NULL, GF = list(ID = NULL, 
    A = NULL), prior = NULL, nIter = 1100, burnIn = 100, thin = 10, 
    thin2 = 1e+10, saveAt = "", minAbsBeta = 1e-09, weights = NULL) 
{

    ETA = NULL
    ETA = list()
    nLT = 0

    cat("This implementation is a simplified interface for the more general\n")
    cat("function BGLR, we keep it for backward compatibility with our package BLR\n")

    warning("thin2 parameter is not used any more and will be deleted in next releases\n",immediate. = TRUE);
   
    cat("Setting parameters for BGLR...\n")
    if (is.null(prior)) {
        cat("===============================================================\n")
        cat("No prior was provided, BGLR will be running with improper priors.\n")
        cat("===============================================================\n")
        prior = list(varE = list(S = NULL, df = 1), varBR = list(S = 0, 
            df = 0), varU = list(S = 0, df = 0), lambda = list(shape = 0, 
            rate = 0, type = "random", value = 50))
    }
    if (!is.null(XF)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XF, model = "FIXED")
    }
    if (!is.null(XR)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XR, model = "BRR", df0 = prior$varBR$df, 
            S0 = prior$varBR$S)
    }
    if (!is.null(XL)) {
        nLT = nLT + 1
        if (prior$lambda$type == "random") {
            if (is.null(prior$lambda$rate)) {
                cat("Setting prior for lambda^2 to beta\n")
                prior$lambda$type = "beta"
                prior$lambda$shape = NULL
                prior$lambda$rate = NULL
            }
            else {
                cat("Setting prior for lambda^2 to gamma\n")
                prior$lambda$type = "gamma"
                prior$lambda$max = NULL
                prior$lambda$shape1 = NULL
                prior$lambda$shape2 = NULL
            }
        }
        ETA[[nLT]] = list(X = XL, model = "BL", type = prior$lambda$type, 
                          rate = prior$lambda$rate, shape = prior$lambda$shape, 
                          max = prior$lambda$max, shape1 = prior$lambda$shape1, 
                          shape2 = prior$lambda$shape2, lambda = prior$lambda$value,
                          minAbsBeta=minAbsBeta)
    }

    #NOTE: In original BLR IDS are used to buid A matrix Z,
    #and then run the model y=Zu+e, u~MN(0,varU*A), and then using the Cholesky factorization
    #it was possible to fit the model. The algorithm used here is different (Orthogonal variables)
    #and it may be that the IDS are not longer necessary

    if (!is.null(GF[[1]])) {
        nLT = nLT + 1
        ETA[[nLT]] = list(K = GF$A, model = "RKHS", df0 = prior$varU$df, 
            S0 = prior$varU$S)
        warning("IDs are not used any more and will be deleted in next releases...\n",immediate. = TRUE) 
    }

    cat("Finish setting parameters for BGLR\n")
    cat("Fitting model using BGLR...\n")
    out = BGLR(y = y, ETA = ETA, df0 = prior$varE$df, S0 = prior$varE$S, 
               nIter = nIter, burnIn = burnIn, thin = thin, saveAt = saveAt, 
               weights = weights)

    #Backward compatibility with BLR
    if (nLT > 0) {
        for (j in 1:nLT) {
            if (ETA[[j]]$model == "FIXED") {
                out$bF = out$ETA[[j]]$b
                out$SD.bF = out$ETA[[j]]$SD.b
            }
            if (ETA[[j]]$model == "BL") {
                out$bL = out$ETA[[j]]$b
                out$SD.bL = out$ETA[[j]]$SD.b
                out$lambda = out$ETA[[j]]$lambda
            }
            if (ETA[[j]]$model == "BRR") {
                out$bR = out$ETA[[j]]$b
                out$SD.bR = out$ETA[[j]]$SD.b
                out$varBR = out$ETA[[j]]$varB
                out$SD.bR = out$ETA[[j]]$SD.varB
            }
            if (ETA[[j]]$model == "RKHS") {
                out$u = out$ETA[[j]]$u
                out$SD.u =out$ETA[[j]]$SD.u
                out$varU = out$ETA[[j]]$varU
            }
        }
    }
    out$ETA = NULL
    class(out) = "BLR"
    return(out)
}

#---------------------------------------------------------------------END OF BGLR FUNCTION-------------------------------------------------------------------------------------#


#Step 8. 
#Running the BGLR with a loop. The function will read the files from the working directory.
#
#Function: BglrLoop()
#
#Value: This is just a loop function for analyzing separate chromosomes. An R .rda save file(s) that
#contain the results of the analysis in a list file object fm and can be 
#loaded into R with load() function.
#
#
#x						Name of the genotype file minus the chromosome number. The genotype file should be 
#						in the format created in step 5.
#
#y						The file extension for the genotype file. The x and y parameters should form the file name
#						by the formula: "x" + chromosome number + "y" i.e. "chromosome1.raw"
#
#out					Desired output file name or directory path and file name
#
#phen					Phenotype to be analysed as vector. 
#
#model					Model to be used in analyses e.g. BayesA, see BGLR package for all available options.				
#
#chrnumbs				vector of the chromosomes to be separated from the data. (for example chrnumbs=c(1:29))
#
#nIter					number of iterations
#
#burnIn					number of cycles burned
#
#Modified				Will the modified version of the BGLR() function be run if TRUE then you need to define the last three arguments as well.
#
#map					If Modified=TRUE, this argument defines the map object for the modified BGLR() function.
#
#CustomLimit			If Modified=TRUE, this argument defines the custom limit for the  modified BGLR() function.
#
#Windowed				Will a windowed analysis also be performed. The windows will be very short, or only the adjacent markers.



BglrLoop<-function(x, Y, out, phen, model, chrnumbs, nIter, burnIn, Modified=FALSE, map=NULL, CustomLimit=4, Windowed=FALSE){

chromnumbers<-chrnumbs

repeat {

kromosomi<-paste(x,chromnumbers[1],Y, sep="")

genot<-read.table(kromosomi, sep="\t", row.names=1, header=TRUE)
save<-paste(out, chromnumbers[1], ".rda", sep="")

n<-nrow(genot)
p<-ncol(genot)
X<-as.matrix(genot)

X<-scale(X,scale=TRUE,center=TRUE)

y<-phen  
y<-scale(y, center=TRUE, scale=TRUE)
y<-as.vector(y)
y<-as.numeric(y)

ETA<-list(MRK=list(X=X,model=model))
if(Modified==FALSE){
	fm<-BGLR(y=y, ETA=ETA, nIter = nIter, burnIn = burnIn,)
	} else {
	
	fm<-BGLR(y=y, ETA=ETA, nIter = nIter, burnIn = burnIn, BigMatrixName=NULL, map=map, CustomLimit=4, Windowed=TRUE)

}
save(fm, file=save)



chromnumbers<-chromnumbers[-1]
if(length(chromnumbers)==0)

break

}

}

#Example
#>library(BGLR)
#>BglrLoop(x="chromosome", Y=".raw", out="BayesC_Phenotype1", phen=phen, model="BayesC", chrnumbs=c(1:10), nIter=100000, burnIn=5000)
#


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#Other Functions

#Producing easily manageable result summary files.
#
#Function: summarize()
#Value: A tab delimited summary file containing information about the analysis. First column will have the chromosome number, 
#second will contain the marker names, the third has absolute effect for that marker and fifth the effect squared.
#Sixth column has the base pair positions for the markers and the seventh column an generic numeric value for each marker
#to make sorting operations easier to manage. 
#
#map			Map file object for the whole genome or, in other words the map file also used in step 1.			
#
#filename		Name of the file(s) to be summarized minus the chromosome number and file extension (.rda)
#				Can be missing if generic=TRUE.
#
#chrnumbs		Vector of the chromosomes to be summarized. (for example chrnumbs=c(1:29)).
#
#generic		If generic=TRUE the function will scan the working directory for files ending in .rda and then 
#				will read and summarize them in reverse chromosomal order according to the vector in chrnumb.
#
#BorC			If model is either Bayes B or Bayes C value of this can be true. The function will then print also the model 
#				selection indicator values that can be used to control Bayes FDR.


summarize<-function(map, filename, chrnumbs, generic=FALSE, BorC=FALSE) {

	map<-map
	
	if(generic==FALSE){
	
	repeat {
		
		writtenin<-paste(filename, chrnumbs[1], ".summary", sep="")
		filename<-paste(filename, chrnumbs[1], ".rda", sep="")

		load(filename)

		gHat<-fm$ETA$MRK$b
		gHat2<-gHat^2
		markers<-as.vector(names(gHat))
		markers<-gsub(".", "-", markers, fixed=TRUE)
		positions<-map[map[,2] %in% noquote(markers),]
		genericvec<-(1:nrow(positions))
		results<-cbind(positions[,1:2], gHat, gHat2, positions[,4], genericvec)
		if(BorC==TRUE){
			modlind<-as.vector(fm$ETA$MRK$d)
			results<-cbind(results,modlind)
			colnames(results)<-c("Chromosome", "Marker", "Effect", "Effect^2", "Position(bp)", "Generic", "ModlSelctInd")
			} else {
			colnames(results)<-c("Chromosome", "Marker", "Effect", "Effect^2", "Position(bp)", "Generic")
		}
		write.table(results, file=writtenin, sep="\t", quote=FALSE, row.names=FALSE)
		
		chrnumbs<-chrnumbs[-1]

		if(length(chrnumbs)==0)

	break

		}
	}
	
	if(generic==TRUE){
		
		files<-list.files(pattern=".rda", full.names=TRUE)
		files2<-list.files(pattern=".rda")
		
		repeat {
		
		filename<-files[grep(paste(tail(chrnumbs, n=1), ".rda", sep=""), files)]
		filename2<-files2[grep(paste(tail(chrnumbs, n=1), ".rda", sep=""), files2)]
		writtenin<-gsub(".rda", ".summary", filename2, fixed=TRUE)
		
		load(filename)
		
		gHat<-fm$ETA$MRK$b
		gHat2<-gHat^2
		markers<-as.vector(names(gHat))
		markers<-gsub(".", "-", markers, fixed=TRUE)
		positions<-map[map[,2] %in% noquote(markers),]
		genericvec<-(1:nrow(positions))
		results<-cbind(positions[,1:2], gHat, gHat2, positions[,4], genericvec)
		if(BorC==TRUE){
			modlind<-as.vector(fm$ETA$MRK$d)
			results<-cbind(results,modlind)
			colnames(results)<-c("Chromosome", "Marker", "Effect", "Effect^2", "Position(bp)", "Generic", "ModlSelctInd")
			} else {
			colnames(results)<-c("Chromosome", "Marker", "Effect", "Effect^2", "Position(bp)", "Generic")
		}
		write.table(results, file=writtenin, sep="\t", quote=FALSE, row.names=FALSE)
		
		chrnumbs<-chrnumbs[-tail(chrnumbs, n=1)]
		files<-files[!files %in% filename]
		files2<-files2[!files2 %in% filename2]
		if(length(chrnumbs)==0)
		
		break
		}	
	}
}




#BFDR
#
#function: BFDR()
#
#value: BFDR -test for the results in files created by summarize() function.
#
#filename				The results file made by summarize() -function as an object.
#
#chrnumbs				This argument is only needed if some of the .summary files are to be
#						tested. If you only want to test some of them and summary files are in the
#						form of ...N.summary e.g. chromosome6.summary you can make this argument
#						a vector of the chromosome numbers for the chromosomes this function is 
#						applied to.
#
#generic				If this argument is true, the working directory will be scanned for ALL
#						.summary files and all of them will be tested.
#
#write					Will the results be printed in the working directory?
#
#z						Significance limit for the test.
#
#


BFDR<-function(filename, chrnumbs=FALSE, generic=TRUE, write=TRUE, z=0.05){
	
	result=0
	IO<-1
	
	if(generic==FALSE){
	
	repeat {
		
		writtenin<-paste(filename, chrnumbs[1], ".result")
		filename<-paste(filename, chrnumbs[1], ".summary", sep="")
		test<-read.table(filename, sep="\t", header=TRUE)
		test<-test[order(test[,7], decreasing=TRUE),]
		test[,7]<-1-test[,7]
	
		count=0
		
		
		for(i in 1:nrow(test)) {
	
		count<-sum(test[1:i,7])/i
			if(count>z){
				n=i-1
				cat(" There were", n, "significant markers in chromosome", chrnumbs[1], "\n")
				test<-test[1:n,]
				break
			}
		}
		
		if(IO==1){
			result<-test
			IO<-0
		} else {
			result<-rbind(result,test)
		}
		chrnumbs<-chrnumbs[-1]
		if(length(chrnumbs)==0)
		break
	}
	}
	
	if(generic==TRUE){
		
	files<-list.files(pattern=".summary", full.names=TRUE)
	files2<-list.files(pattern=".summary")
	d=1
	
	repeat {
		
		filename<-files[grep(paste(tail(chrnumbs, n=1), ".summary", sep=""), files)]
		
		if(d==1){
		filename2<-files2[grep(paste(tail(chrnumbs, n=1), ".summary", sep=""), files2)]
		filename2<-gsub("[[:digit:]]+", "", filename2)
		writtenin<-gsub(".summary", ".result", filename2, fixed=TRUE)
		d=0
		}
	
		test<-read.table(filename, sep="\t", header=TRUE)
		test<-test[order(test[,7], decreasing=TRUE),]
		test[,7]<-1-test[,7]
	
		count=0
		
		for(i in 1:nrow(test)) {
	
		count<-sum(test[1:i,7])/i
			if(count>z){
				n=i-1
				cat(" There were", n, "significant markers in chromosome", tail(chrnumbs, n=1), "\n")
				test<-test[1:n,]
				break
			}
		}
		test<-test[order(test$Position.bp.),]
		if(IO==1){
			result<-test
			IO<-0
		} else {
			result<-rbind(result,test)
		}
		chrnumbs<-chrnumbs[-tail(chrnumbs, n=1)]
		files<-files[!files %in% filename]
		if(length(chrnumbs)==0)
		
		break
	}
	
	Generic2<-c(1:nrow(result))
	result<-result[order(result$Chromosome),]
	result<-cbind(result,Generic2)
	
	if(write==TRUE){
		write.table(result, file=writtenin, quote=FALSE, sep="\t", row.names=FALSE)
		cat(" The result file: \n", writtenin, "\n", "is in the working directory: \n", paste(getwd()), "\n")
	} else {
		return(result)
	}
	}	
}







	
