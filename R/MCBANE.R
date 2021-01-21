
#the codetools package provides the checkUsage function which ensures that the given function is self contained
library(codetools)


#I don't think this was used. I believe it was similar to the analysis below, but used ANOVA instead of T-tests to look for differences in expression between alleles
testAllSNPsANOVA = function(tagData, sample, allele="allele", snpid="snpid", tag="tag",probeid="probeid", primarySNP="primarySNP", logFC="logFC", minTags=5){
  message("Filtering innapropriate SNPs for ANOVA test");
  tagData = tagData[!is.na(tagData[[logFC]]),] # remove tags without a logFC (e.g. bad data)
  allelesPerSNP = cast(unique(tagData[c(allele,snpid)]), as.formula(sprintf("%s ~ .",snpid)), fun.aggregate = function(x){length(unique(x))}, value=allele)
  names(allelesPerSNP)[ncol(allelesPerSNP)]="numAlleles";
  allelesPerSNP = allelesPerSNP[allelesPerSNP$numAlleles>2,]
  nonPrimarySNPs = cast(unique(tagData[c(primarySNP,snpid)]), as.formula(sprintf("%s ~ .",snpid)), fun.aggregate = function(x){length(unique(x))}, value=primarySNP)
  
  names(nonPrimarySNPs)[ncol(nonPrimarySNPs)]="nonPrimary";
  nonPrimarySNPs = nonPrimarySNPs[nonPrimarySNPs$nonPrimary>1,]
  eitherSNPs = unique(c(as.character(allelesPerSNP[[snpid]]), as.character(nonPrimarySNPs[[snpid]])))
  message(sprintf("Running ANOVA for %i SNPs",length(eitherSNPs)))
  tagData = tagData[tagData[[snpid]] %in% eitherSNPs,]                 
  tagData$distToPrimary = abs(tagData$posPrimary - tagData$position);
  tagData = tagData[tagData$distToPrimary<=60,];
  allSNPs = unique(tagData[[snpid]]);
  allComparisons = data.frame();
  allSamples = unique(tagData[[sample]])
  for (snpi in 1:length(allSNPs)){
    snp = allSNPs[snpi];
    message(sprintf("snp %i/%i: ",snpi,length(allSNPs)), allSNPs[snpi]);
    curData = tagData[tagData[[snpid]] == snp,];
    for(sampleID in allSamples){
      #message(sampleID)
      curData2 = curData[curData[[sample]]==sampleID,];
      #message(str(curData2))
      if (((!all(curData2[[snpid]]==snp)) || length(unique(curData2[[allele]]))>2)){
        #there are multple SNPs present on this probe; a haploblock || there are more than two alleles:
        #do (M)ANOVA test
        curProbeTagVals = unique(curData2[c("probeid","tag","logFC")]);
        curFactors = cast(unique(curData2[c("probeid","snpid","primarySNP","allele")]), formula = probeid +primarySNP ~ snpid, value = "allele")
        if (sum(is.na(curFactors))>0){  # replace NAs with a new factor level: "missing"
          for ( colWithNAs in (1:ncol(curFactors))[apply(is.na(curFactors),2,any)]){
            levels(curFactors[[colWithNAs]]) = c(levels(curFactors[[colWithNAs]]),"missing");
            curFactors[is.na(curFactors[[colWithNAs]]),colWithNAs] = "missing";
          }
        }
        curFactors = as.data.frame(curFactors);
        #unique(curFactors[!apply(is.na(curFactors),2,any)]) # remove any with NAs since these are confounded with primarySNP
        curFactors=unique(curFactors[apply(curFactors,2,function(x){length(unique(x))})>1])  ## remove any without levels
        if (snp %in% colnames(curFactors)){ # continue only if there was more than one allele for the SNP tested
          curSNPs = names(curFactors)[2:ncol(curFactors)];
          names(curFactors)[2:ncol(curFactors)] = gsub(":","_",paste("x",names(curFactors)[2:ncol(curFactors)], sep=""))
          curSNPAliases=names(curFactors)[2:ncol(curFactors)];
          curAOVData = merge(curProbeTagVals, curFactors, by = c("probeid"));
          
          if (nrow(curAOVData) >= 10){
            curAOV = aov(as.formula(sprintf("logFC ~ %s", paste(curSNPAliases, collapse=" + "))),  data = curAOVData)
            tukeyRes = TukeyHSD(curAOV, which = curSNPAliases[curSNPs==snp])
            if (!is.na(names(tukeyRes))){
              for (k in 1:nrow(tukeyRes[[1]])){
                alA = gsub("(.*)-(.*)","\\1",row.names(tukeyRes[[1]])[k])
                alB = gsub("(.*)-(.*)","\\2",row.names(tukeyRes[[1]])[k])
                allComparisons = rbind(allComparisons, data.frame(statType = "ANOVA", P = tukeyRes[[1]][k,4], diff = tukeyRes[[1]][k,1], alleleA = alA, alleleB = alB, snp = snp, nA = sum(curData2[[snpid]]==snp & curData2[[allele]]== alA), nB =  sum(curData2[[snpid]]==snp & curData2[[allele]]== alB), sample = sampleID));
              }
            }
          }
        }
      }
    }
  }
  return(allComparisons)
}

checkUsage(testAllSNPsANOVA)

#' This is the main function that calculates ASE by either t-test or rank sum test
#' @param tagData data.frame containing the tag-level data, including (at a minimum) log FC and the fields referred to in the parameters of this function
#' @param sample the name of the column containing the sample IDs for each experiment
#' @param allele the name of the column containing the allele that each tag is parcoding
#' @param snpid the name of the column containing the SNP ID (e.g. rs58792)
#' @param tag the name of the column containing the tag sequence
#' @param probeid the name of the column containing the probe ID
#' @param logFC the name of the column containing the log(RNA/DNA) values for each tag
#' @param test the type of test. "t.test" for Student's t-test; "ranksum" for the Wilcoxon rank sum test.
#' @param minTags the min number of tags per SNP that are required before a test is performed
#' @return Returns a data.frame containing the ASE results
testAllSNPs = function(tagData, sample, allele="allele", snpid="snpid", tag="tag",probeid="probeid", logFC="logFC", test="t.test", minTags=5){
  if (!(test %in% c("t.test","ranksum"))){
    stop("Unrecognized value for 'test'.  Must be one of c('t.test','ranksum')")
  }
  tagData = tagData[!is.na(tagData[[logFC]]),] # remove tags without a logFC (e.g. bad data)
  #TODO: we used to do this, but maybe this would be better to do somewhere else. We might not have this info and why did we include the SNPs if we didn't plan to study them?
  #tagData$distToPrimary = abs(tagData$posPrimary - tagData$position);
  #tagData = tagData[tagData$distToPrimary<=60,];
  allSNPs = unique(tagData[[snpid]]);
  allComparisons = data.frame();
  allSamples = unique(tagData[[sample]])
  for (snpi in 1:length(allSNPs)){
    snp = allSNPs[snpi];
    message(sprintf("snp %i/%i: ",snpi,length(allSNPs)), allSNPs[snpi]);
    curData = tagData[tagData[[snpid]] == snp,];
    for(sampleID in allSamples){
      #message(sampleID)
      curData2 = curData[curData[[sample]]==sampleID,];
      #message(str(curData2))
      #do t.test or rank sum
      if (test %in% c("t.test", "ranksum")){
        curAlleles = unique(curData2[[allele]])
        #message(paste(curAlleles, collapse=", "))
        if (length(curAlleles)>=2){
          for (a1i in 1:(length(curAlleles)-1)){ # compare all pairs of alleles
            for(a2i in (a1i+1):length(curAlleles)){
              if (sum(curData2[[allele]]==curAlleles[a1i])>=minTags && sum(curData2[[allele]]==curAlleles[a2i]) >=minTags){ # min number of tags to bother comparing
                if (test=="t.test"){
                  #message(sprintf("a1 = %i, a2= %i", sum(curData2[[allele]]==curAlleles[a1i]), sum(curData2[[allele]]==curAlleles[a2i])))
                  curTT = t.test(curData2[[logFC]][curData2[[allele]]==curAlleles[a1i]], curData2[[logFC]][curData2[[allele]]==curAlleles[a2i]]);
                  allComparisons = rbind(allComparisons, data.frame(statType = "t.test", P = curTT$p.value, diff=curTT$estimate[1]-curTT$estimate[2], alleleA = curAlleles[a1i],alleleB = curAlleles[a2i], snp = snp, nA=sum(curData2[[allele]]==curAlleles[a1i]), nB = sum(curData2[[allele]]==curAlleles[a2i]), meanA=mean(curData2[[logFC]][curData2[[allele]]==curAlleles[a1i]]), meanB=mean(curData2[[logFC]][curData2[[allele]]==curAlleles[a2i]]), sample = sampleID));
                }else{
                  curRS = wilcox.test(curData2[[logFC]][curData2[[allele]]==curAlleles[a1i]], curData2[[logFC]][curData2[[allele]]==curAlleles[a2i]]);
                  allComparisons = rbind(allComparisons, data.frame(statType = "ranksum", P = curRS$p.value, diff=median(curData2[[logFC]][curData2[[allele]]==curAlleles[a1i]]) - median(curData2[[logFC]][curData2[[allele]]==curAlleles[a2i]]), alleleA = curAlleles[a1i],alleleB = curAlleles[a2i], snp = snp, nA=sum(curData2[[allele]]==curAlleles[a1i]), nB = sum(curData2[[allele]]==curAlleles[a2i]), medianA=median(curData2[[logFC]][curData2[[allele]]==curAlleles[a1i]]), medianB=median(curData2[[logFC]][curData2[[allele]]==curAlleles[a2i]]), sample = sampleID));
                }
              }
            }
          }
        }
      }
    }
  }
  return(allComparisons)
}
checkUsage(testAllSNPs)


#' This is the function for inputting the MPRA data. Here, each file contains headers but no row names. 
#' @param IDs a vector of sample IDs, each corresponding to the entries in files
#' @param files a vector of file paths, one for each ID in IDs
#' @return Returns a list where each ID points to a vector of counts.
readMPRAData = function(IDs, files){
  allData = data.frame();
  i=1;
  for(i in 1:length(IDs)){
    curData = read.table(files[i],stringsAsFactors=F, row.names=NULL, header=T);
    if (i==1){
      names(curData)[3]=IDs[i];
      allData = curData;
    }else{
      allData[IDs[i]] = curData$count;
    }
  }
  return(allData)
}
checkUsage(readMPRAData)

#' This is the function for normalizing the GC bias present with MPRA tags. Tag expression correlated with GC content and this corrects for this by calculating the mean expression for each GC bin, and subtracting this from tags within that bin.
#' @param logFCMat a matrix of tag log(RNA/DNA) values. Columns are samples, rows are tags.
#' @param tags the tag sequences (detaults to row.names(logFCMat))
#' @param makePlot Logical. Whether or not to return a ggplot plot p which will display the normalization.
#' @return Returns a list containing the normalized logFC matrix, the normalization factors themselves, and (if makePlot is TRUE), the normalization plot.
normTagsByGC=function(logFCMat, tags = row.names(logFCMat), makePlot=T){
  message("Counting GC content")
  gcContent = str_count(tags, "[GC]");
  uGC = sort(unique(gcContent));
  #normalize the data by subtracting the mean from each GC bin - do it all samples together
  message("Calculating normalization factors");
  curNormFactors = data.frame(gc = uGC, norm =rep(0,length(uGC)), n=NA);
  for(g in 1:length(uGC)){
    curNormFactors$norm[g] = mean(logFCMat[gcContent==uGC[g],], na.rm = T);
    curNormFactors$n[g] = sum(!is.na(logFCMat[gcContent==uGC[g],]))
  }
  if (makePlot){
    message("Plotting data")
    #print pre-norm for one example
    logFCMat2 = as.data.frame(logFCMat)
    logFCMat2$gc = gcContent; logFCMat2$tag=tags;
    logFCMat2 = as.data.frame(melt(logFCMat2, id.vars=c("gc","tag")));
    names(logFCMat2)[(ncol(logFCMat2)-1):ncol(logFCMat2)]=c("sample","logFC");
    logFCMatSummary = cast(logFCMat2, gc+sample ~ ., value="logFC", fun.aggregate = function(x){list(mean=mean(x, na.rm=T), sd = sd(x,na.rm=T))})
    #return(logFCMat2)
    p = ggplot(logFCMatSummary, aes(y=mean, x= gc)) + theme_classic()+geom_point()+facet_grid(sample ~.) + geom_point(data=curNormFactors,aes(y=norm),colour="red") + scale_y_continuous(expand=c(0,0))+geom_hline(yintercept = 0); print(p);
  }
  #subtract mean in each bin
  message("normalizing");
  for(g in 1:length(uGC)){
    logFCMat[gcContent==uGC[g],] = logFCMat[gcContent==uGC[g],] - curNormFactors$norm[g];
  }
  #logFCMat = logFCMat - median(logFCMat, na.rm=T) # unclear why I need this
  if (makePlot){
    return(list(normLogFC = logFCMat, normFactors = curNormFactors, plot=p))
  }else{
    return(list(normLogFC = logFCMat, normFactors = curNormFactors))
  }
}
checkUsage(normTagsByGC)


#' Not sure this function is used.
#' This is the function for calculating the log of the tag data, with the option for including a pseudo count.
#' @param tagData a data.frame of tags (rows) by samples (columns), where entries contain read counts for that tag. Not sure what first two columns are, but these are not counts (e.g. could be tag sequence and probe ID)
#' @param pseudocount a count to add to each read count to avoid log(0). Defaults to 0.5
#' @return Returns the same data frame as tagData, but counts are now log(count+pseudocount)
logTagData=function(tagData, pseudocount=0.5){
  tagData[3:(ncol(tagData))] = log2(pseudocount+tagData[3:(ncol(tagData))]);
  return(tagData)
}
checkUsage(logTagData)


#' This is the function for calculating tag log(RNA/DNA) values for an MPRA experiment.
#' @param tagCounts a data.frame of tag counts, one per sample, and one containing the input tags. Usually input is amplicon sequencing of the MPRA plasmid library. One column is named 'tag' one 'enhancer' (we should update this).
#' @param inputColName the name of the column in tagCounts that contains the input (DNA) tag counts.
#' @param center function with which to centre the data for each sample (default=median)
#' @param minInput the minimum reads required within the input (DNA) before we include this tag (default = 30)
#' @param minOutput the minimum reads required in the output (RNA) before we include this tag (default = 4)
#' @return Returns the logFC as a matrix
getTagFCs=function(tagCounts, inputColName,  center=median, minInput=30, minOutput=4){
  outCountMat = tagCounts[, colnames(tagCounts)!=inputColName];
  outCountMat$enhancer=NULL; outCountMat$tag=NULL
  outCountMat = as.matrix(outCountMat)
  inputCounts = tagCounts[[inputColName]];
  outTotals = apply(outCountMat, 2, sum, na.rm=T)
  inputTotal = sum(inputCounts, na.rm=T);
  badOutputs = outCountMat < minOutput
  badInputs = inputCounts < minInput
  outCountMat = scale(outCountMat, center=F, scale=outTotals/1E6) #to CPM
  inputCounts = inputCounts/inputTotal*1E6 #to CPM
  logFC = outCountMat / inputCounts;
  logFC = log2(logFC);
  #filter out datapoints of insufficient coverage
  message(sprintf("Removing %i tags for poor coverage in output (reads<%i); %f per sample", sum(badOutputs), minOutput, sum(badOutputs)/ncol(outCountMat)))
  message(sprintf("Removing %i tags for poor coverage in input (reads<%i);",sum(badInputs),minInput))
  logFC[badOutputs] = NA;
  logFC[badInputs,] = NA
  row.names(logFC) = tagCounts$tag;
  #center the data, by default using the median
  if (!is.null(center)){
    allSamples = colnames(logFC);
    for ( j in 1:length(allSamples)){ 
      #Theoretically, I would divide each logFC by the ratio of the total counts observed for that sample: total for the input 
      #this is equivalent to: going to subtract (log2(total_RNA) - log2(total_input)) since we already log2 transformed the data
      #print((log2(mean(outCountMat[,colnames(outCountMat)==allSamples[j]])) - log2(mean(inputPreViralCounts))))
      #logFC[,colnames(logFC)==allSamples[j]] = logFC[,colnames(logFC)==allSamples[j]] - (log2(mean(outCountMat[,colnames(outCountMat)==allSamples[j]])) - log2(mean(inputPreViralCounts)));
      #instead, I will subtract the median logFC so that the distributions are centered
      message(sprintf("Scaling %s logFC by %g", allSamples[j], center(logFC[,colnames(logFC)==allSamples[j]], na.rm=T)))
      logFC[,colnames(logFC)==allSamples[j]] = logFC[,colnames(logFC)==allSamples[j]] - center(logFC[,colnames(logFC)==allSamples[j]], na.rm=T)
    }
  }
  return(logFC)
}
checkUsage(getTagFCs)


#' This is the function for determining which tags contain any of a set of given k-mers
#' @param kmers a vector of k-mers
#' @param tags vector of tag sequences
#' @return Returns a logical for whether each tag contains any of the included k-mers
tagContainsKmer = function(kmers, tags){
  containsKmer = rep(F, length(tags));
  for(kmer in kmers){
    containsKmer = containsKmer | grepl(kmer, tags)
  }
  return(containsKmer);
}
checkUsage(tagContainsKmer)



#' This is the function for identifying likely outlier k-mers. These are k-mers which may affect the expression of the tag if present within it. These can include splicing motifs, etc.
#' @param logFCMat Matrix of log(RNA/DNA) values.
#' @param tags tag sequences (default=row.names(logFCMat))
#' #param k The length of k-mers to consider (run time increases with 4^k). (default=5)
#' @return Returns a data.frame containing mean(logFC), sd(logFC), and n (number of probes containing the k-mer) for each k-mer. 
findOutlierKmers = function(logFCMat, tags = row.names(logFCMat), k=5){
  alphabet = c("A","T","G","C")
  kmers = ""
  for(i in 1:k){ #modified from imminfo/tcr
    kmers =  unlist(lapply(kmers, function (kmer) { paste0(kmer, alphabet) }))
  }
  kmers = data.frame(kmer=kmers, mean=NA, sd=NA, n=NA)
  for(i in 1:nrow(kmers)){
    message(sprintf("%i/%i",i,nrow(kmers)))
    subset = grepl(kmers$kmer[i], tags)
    kmers$mean[i] = mean(logFCMat[subset,], na.rm=T)
    kmers$sd[i] = sd(logFCMat[subset,], na.rm=T)
    kmers$n[i] = sum(!is.na(logFCMat[subset,]))
  }
  kmers = kmers[order(kmers$mean),]
  return(kmers)
}
checkUsage(findOutlierKmers)

#' This is the function for removing tags that are NA for one reason or another.
#' @param tagData data.frame (tags and probes in the first two columns) containing read counts for each sample (columns 3+)
#' @return Returns the tagData, but after removing tags with 0 reads in all samples. 
filterMissingTags=function(tagData){
  totalTags = nrow(tagData);
  #throw out any tags that were never observed in any data
  tagData = tagData[apply(tagData[,3:ncol(tagData)]>0,1,any),];
  observedTags = nrow(tagData);
  message(sprintf("%i of %i tags observed (%g%%)",observedTags,totalTags, 100*observedTags/totalTags));
  return(tagData)
}
checkUsage(filterMissingTags)
