       ###############################
       ######## Assignment 10 ########
       ###############################
       
      #Discrepancy Analysis

       rm(list = ls())
       #setwd("/Users/emanuelastruffolino/Desktop/SequenceCourse/Assignments")
       #getwd()
       library (TraMineRextras)
       
       #1. Consider the state sequence object biofam.seq and the matrix dOM of 
            #pairwise OM dissimilarities based on the properties matrix considered 
            #in the previous assignments.
       data(biofam)
       mycol<-brewer.pal(,"RdBu")
       biofam.lab<-c("Parent", "Left", "Married","Left+Marr","Child","Left+Child","Left+Marr+Child","Divorced")
       biofam.shortlab<-c("P", "L", "M", "LM", "C","LC","LMC", "D")
       biofam.seq <- seqdef(biofam[,10:25],states=biofam.shortlab,labels=biofam.lab)
       summary (biofamseq)
       
       properties <- matrix(c(# left, married, child, divorced
         0, 0, 0, 0,  # parent
         1, 0, 0, 0,  # left
         0, 1, .5, 0, # marr
         1, 1, 0, 0,  # left+marr
         0, 0, 1, 0,  # child
         1, 0, 1, 0,  # left+child
         1, 1, 1, 0,  # left+marr+child
         .5, 1, .5, 1 # divorced
       ), 8, 4, byrow=TRUE)
       scm <- as.matrix(dist(properties))
       indel <- .5*max(scm)
       distOM <- seqdist(biofam.seq, method="OM", indel=indel, sm=scm, full.matrix = FALSE)
       
       #2. Create a cohort covariate, separating the individual born after the 
           #second World War from those born before.
       biocoho <- cut(biofam$birthyr, c(1900,1946,1960), labels=c("Before W", "After W"), right=FALSE)
       
       #3. We are interested in studying the rise of new kind of family trajectories. 
          #Plot the sequences according to the cohort factor.
       seqplot(biofam.seq,biocoho)
       seqIplot(biofam.seq,biocoho, sortv="from.end")
       
       #4. To test if the differences are significative, compute the association between 
          #the sequences and the cohort covariate using dissassoc.
       library(WeightedCluster)
       weight<-attr(biofam.seq, "weight")
       coho.assoc <- dissassoc(distOM, group = biocoho, R = 5000,
                                 weights = weight, weight.permutation = "diss")
       print(coho.assoc)
       
       #5. Some have argued that one of the main change is related with the diversification 
          #of the family trajectories (this assumption is called destandardisation of the life course). 
          #To test this hypothesis, we can look at the difference of discrepancies over cohorts. 
          #What do you think? Are the differences significant?
       
       #6. Explore the evolution of the association over time using seqdiff.
       sm <- as.matrix(dist(properties))
       indel <- .5*max(sm)       
       biofam.diff<-seqdiff(biofam.seq,group=biocoho,seqdist_arg = list(method = "OM", indel = 1.5, sm=sm))
       plot( biofam.diff, stat=c("Pseudo R2", "Levene"))
       plot( biofam.diff, stat="discrepancy")
       
       #7. Build a sequences regression tree using the sex, cohort and plingu02 (langue of interview) 
       biofam <- droplevels(biofam)
       tree <- seqtree(biofam.seq ~ sex + biocoho + plingu02,
                        data = biofam, R = 5000, diss = distOM, weight.permutation = "diss")
       
       # plot it using seqtreedisplay.
       seqtreedisplay(tree, type = "r", showtree = TRUE,
                      showdepth = TRUE, dist.matrix=distOM, legend.fontsize=2)

       #this last command doesn't work because: "Error in dist.matrix[ind, ind] : wrong number of dimensions"

