#PURPOSE: Modifications to some of the permulation scripts for murine sperm project
	#Modified getPermsBinary SSM version to set root for each gene tree, because the master tree root species is not always present in every gene tree
	#Faster SSM permulations that do not check for permulated trees that match the true tree structure
	#Quiet versions of functions that do not print out "Species not present..." message every time
	#Incorporated simBinPhenoRank into main wrapper functions - NOT TESTED (I abandoned this because it didn't work w/ SSM perms)
	#Added a function to make a custom "sisters_list" for each gene tree when running SSM permulations
	#Note: library(RERconverge) needs to be loaded for these to work

####getPermsBinaryQuiet()####
#Copied from getPermsBinary() [permulationFunctions.R]
#Modified to set root_sp to first taxa listed in each indivdual gene tree instead of the same root species for every gene tree
#This is to deal with the fact that not every species is in every gene tree due to filtering, missing data, etc
#Elysia says the outgroup/root species is unimportant for permulations, so it is okay to set a different root_sp for each gene tree, and to set the outgroup arbitrarily this way
#Also runs the "quiet" versions of functions that DO NOT print "Species not present..." message

#'Calculates permuted correlation and enrichment statistics for binary phenotype
#' @param numperms An integer number of permulations
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param root_sp The species to root the tree on
#' @param RERmat An RER matrix calculated using \code{\link{getAllResiduals}}.
#' @param trees treesObj from \code{\link{readTrees}}
#' @param mastertree A rooted, fully dichotomous tree derived from the treesObj master tree from \code{\link{readTrees}}.  Must not contain species not in traitvec
#' @param permmode Mode of binary permulation ("cc" for Complete Cases (default), "ssm" for Species Subset Match)
#' @param method statistical method to use for correlations (set to "k" (default) for Kendall Tau test)
#' @param min.pos minimum number of foreground species (default 2)
#' @param trees_list A list containing the trees of all genes of interest (formatted like trees in treesObj from \code{\link{readTrees}})
#' @param calculateenrich A boolean variable indicating if null permulation p-values for enrichment statistics
#' @param annotlist Pathway annotations
#' @return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#' @export
getPermsBinaryQuiet=function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees, mastertree, permmode="cc", method="k", min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL){
  pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels

  if (permmode=="cc" | permmode=="rank"){
    print("Running CC permulation")

    print("Generating permulated trees")
    if (permmode=="cc"){
        print("Running CC permulation")
        permulated.binphens = generatePermulatedBinPhenQuiet(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    } else {
        print("Running CC permulation using simBinPhenoRank")
	permulated.binphens = list()
	for(i in 1:numperms){
	        permulated.binphens[[length(permulated.binphens)+1]] = simBinPhenoRank(trees, root_sp, fg_vec, sisters_list, plotTreeBool=F)
    	}
    }
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])

    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)

    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)

    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }

  } else if (permmode=="ssm"){
    print("Running SSM permulation")

    if (is.null(trees_list)){
      trees_list = trees$trees
    }

    RERmat = RERmat[match(names(trees_list), rownames(RERmat)),]

    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatchedQuiet(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec)

    # Get species membership of the trees
    print("Getting species membership of the trees")
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=mastertree,foregrounds=fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names

    colnames(df.converted) = col_names
    rownames(df.converted) = row_names

    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)

    spec.members = df.converted$spec.members

    # Group gene trees based on the similarity of their species membership
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups

    # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
    unique.trees = trees_list[ind.unique.trees]

    # precompute clade mapping for each unique tree
    print("Computing clade mappings for each unique tree")
    unique.map.list = mapply(matchAllNodesClades,unique.trees,MoreArgs=list(treesObj=trees))

    # calculate paths for each permulation
    print("Calculating paths for each permulation")
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens,unique.map.list,trees)

    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)){
      permulated.paths[[j]] = vector("list",length=numperms)
    }
    for (i in 1:length(unique.permulated.paths)){
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)){
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names

    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i,])
    corMatList = mapply(calculateCorPermuted,permulated.paths,RERmat.list)
    permPvals = extractCorResults(corMatList,numperms,mode="P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList,numperms,mode="Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals)*-log10(permPvals)
    rownames(permStatvals) = names(trees_list)

  } else {
    stop("Invalid binary permulation mode.")
  }

  if (calculateenrich){
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, trees, plotTree=F)
    realpaths = tree2PathsCladesQuiet(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method=method, min.pos=min.pos)
    realstat =sign(realresults$Rho)*-log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, outputGeneVals=F)

    #sort real enrichments
    groups=length(realenrich)
    c=1
    while(c<=groups){
      current=realenrich[[c]]
      realenrich[[c]]=current[order(rownames(current)),]
      c=c+1
    }
    #make matrices to fill
    permenrichP=vector("list", length(realenrich))
    permenrichStat=vector("list", length(realenrich))
    c=1
    while(c<=length(realenrich)){
      newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich[[c]])))
      rownames(newdf)=rownames(realenrich[[c]])
      permenrichP[[c]]=newdf
      permenrichStat[[c]]=newdf
      c=c+1
    }

    counter=1;
    while (counter <= numperms){
      stat = permStatvals[,counter]
      names(stat) = rownames(RERmat)
      enrich=fastwilcoxGMTall(na.omit(stat), annotlist, outputGeneVals=F)
      #sort and store enrichment results
      groups=length(enrich)
      c=1
      while(c<=groups){
        current=enrich[[c]]
        enrich[[c]]=current[order(rownames(current)),]
        enrich[[c]]=enrich[[c]][match(rownames(permenrichP[[c]]), rownames(enrich[[c]])),]
        permenrichP[[c]][,counter]=enrich[[c]]$pval
        permenrichStat[[c]][,counter]=enrich[[c]]$stat
        c=c+1
      }
      counter = counter+1
    }
  }

  if(calculateenrich){
    data=vector("list", 5)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    data[[4]]=permenrichP
    data[[5]]=permenrichStat
    names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  } else {
    data=vector("list", 3)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    names(data)=c("corP", "corRho", "corStat")
  }
  data
}

#NOTE: Below are some functions I didn't change but copied because they are marked as "internal" so the functions that call them can't find them unless they're in the same file

#' @keywords internal
getForegroundsFromBinaryTree=function(tree){
  nameEdgesPerms.tree = nameEdgesPerms(tree)
  edge.length = as.logical(tree$edge.length)
  foregrounds = nameEdgesPerms.tree[edge.length]
  ind.tip = which(foregrounds != "")
  foregrounds = foregrounds[ind.tip]
  return(foregrounds)
}

#' @keywords internal
nameEdgesPerms=function(tree){
  if (is.null(tree$tip.label)) {
    nn = NULL
  } else {
    nn=character(nrow(tree$edge))
    iim=match(1:length(tree$tip.label), tree$edge[,2])
    nn[iim]=tree$tip.label
  }
  nn
}

#NOTE: I didn't change anything in getForegoundInfoClades; I just copied it here to add print statements to figure out where permulations were stalling
#'Generates a binary phenotype tree and foreground clades information using the list of tip foreground animals, the presence of foreground common ancestors, and their phylogenetic relationships
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param trees treesObj from \code{\link{readTrees}}
#' @param plotTree A boolean indicator for plotting the output tree (default=FALSE)
#' @param useSpecies An array containing the tip labels in the output tree
#' @return output.list A list containing 1) "tree" = a binary phenotype tree corresponding to the input information, 2) "fg.sisters.table" = a table containing all sister species in the foreground set
#' @export
getForegroundInfoClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,useSpecies=NULL){
  if (length(useSpecies)==0){
    useSpecies = trees$masterTree$tip.label
  }

  if (is.null(sisters_list)){
    fg.sisters.table=NULL
    fg_tree = foreground2Tree(fg_vec, trees, plotTree=F, clade="terminal", useSpecies=useSpecies)
  } else {
    # Start with a temp phenotype tree assuming that all internal nodes are foregrounds
    #print("Starting foreground2Tree")
    fg_tree = foreground2Tree(fg_vec,trees,plotTree=F,clade="all",useSpecies=useSpecies)
    #print("Done with foreground2Tree")
    #write.tree(fg_tree, "temp.tre", append=T)
    edge = fg_tree$edge
    edge.length=fg_tree$edge.length

    ind.fg.edge = which(edge.length == 1)
    #print(ind.fg.edge)
    nodeIds.fg.edge = edge[ind.fg.edge,] # all foreground edges in the temp tree
    #print(nodeIds.fg.edge)

    tip.sisters = vector("integer",length=0)
    for (i in 1:length(sisters_list)){
      sisters = sisters_list[[i]]
      #print(sisters)
      nodeId.sisters = which(useSpecies %in% sisters)
      if (length(nodeId.sisters)>0){
        tip.sisters = c(tip.sisters,nodeId.sisters)
      }
    }
    #print(tip.sisters)

    # Find and correct the pairs
    fg.sisters.table = matrix(nrow=0,ncol=2)
    colnames(fg.sisters.table) = c("species1","species2")
    if (length(as.vector(nodeIds.fg.edge)) > 2){
      all.nodeId.ca = sort(nodeIds.fg.edge[,1])
      count_all_nodeId_ca = table(all.nodeId.ca)
      unq.nodeId.ca = unique(all.nodeId.ca)
      fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
      nodes_addressed = NULL
      #print("Starting while loop")
      #print(unq.nodeId.ca)
      while (length(unq.nodeId.ca) != length(nodes_addressed)){
        nodeId.ca = sort(all.nodeId.ca[which(!(all.nodeId.ca %in% nodes_addressed))])
        #print(nodeId.ca)
        if (length(nodeId.ca) == 1){
          nodes_addressed = c(nodes_addressed, nodeId.ca)
        } else {
          for (nn in 1:(length(nodeId.ca)-1)){

            if (nodeId.ca[nn] == nodeId.ca[nn+1]){
              nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
              #print(nodeId.desc)  
            if (length(which(nodeId.desc %in% tip.sisters)) > 0){
                fg_ca = c(fg_ca,nodeId.ca[nn])
                #print(fg_ca)
                fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
                #print(fg.sisters.table)
                #print("here1")
                nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              } else {
                if (length(which(fg_tree$tip.label[nodeId.desc] %in% fg_vec)) == 2){
                  fg_tree$edge.length[which(edge[,2]==nodeId.ca[nn])] = 0
                  #print("here2")
                  nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                } else {
                  if (length(which(nodeId.desc %in% nodes_addressed)) == 2){
                    fg_ca = c(fg_ca,nodeId.ca[nn])
                    fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
                    #print("here3")
                    nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                  }
                  #print("no else here...")
                  #print(length(which(nodeId.desc %in% nodes_addressed)))
                }
              }
            } else {
              nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
              if (length(nodeId.desc) == 2){
                if (nodeId.ca[nn] != nodeId.ca[nn-1]){
                  fg_tree$edge.length[which(edge[,2] == nodeId.ca[nn])] = 0
                  #print("here4") 
                  nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                  nodes_addressed = unique(nodes_addressed)
                }
              } else {
                #print("here5")
                nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              }
            }
          }
        }
        #print(nodes_addressed)
      }
      rownames(fg.sisters.table) = fg_ca
    }
  }
  if (plotTree==T){
    plot(fg_tree)
  }
  output.list = list("fg.sisters.table"=fg.sisters.table,"tree"=fg_tree)
  output.list
}

#'Produces one binary permulation based on ranking of simulated branch lengths
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec a vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @return A binary permulated tree
#' @export
simBinPhenoRankQuiet=function(trees, root_sp, fg_vec, sisters_list=NULL, plotTreeBool=F){
  #print("In simBinPhenoRank")
  mastertree = trees$masterTree
  #print(mastertree)
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
  fg_tree = res$tree
  pathvec = tree2PathsCladesQuiet(fg_tree, trees)

  t=root.phylo(trees$masterTree, root_sp, resolve.root = T)
  ratem=ratematrix(t, pathvec)


  x = rnorm(n=length(t$edge.length))
  sd = sqrt(as.vector(ratem)*t$edge.length)

  y = matrix(0,nrow(t$edge),ncol(t$edge))
  alpha=0.1
  n = length(t$tip)
  for(i in 1:length(x)){
    if(t$edge[i,1]==(n+1))
      y[i,1]<-alpha # if at the root
    else
      y[i,1]<-y[match(t$edge[i,1],t$edge[,2]),2]
    y[i,2]<-y[i,1]+x[i]
  }

  rm(x)

  x<-c(y[1,1],y[,2])
  names(x)<-c(n+1,t$edge[,2])

  x<-x[as.character(1:(n+t$Nnode))]

  simphentree = t
  simedge = x[simphentree$edge[,2]]
  simphentree$edge.length = unname(simedge) + abs(min(simedge)) + 0.001

  numfg= sum(fg_tree$edge.length)

  simedgesort = sort(simphentree$edge.length, decreasing=T)
  simphenthreshold = simedgesort[numfg]

  bmphentree = simphentree
  simbinedge = simphentree$edge.length
  simbinedge[which(simbinedge < simphenthreshold)] = 0
  simbinedge[which(simbinedge >= simphenthreshold)] = 1
  bmphentree$edge.length = simbinedge

  if(plotTreeBool){
    plot(bmphentree)
  }

  return(bmphentree)
}

####generatePermulatedBinPhenQuiet()####
#EK - Copied from generatePermulatedBinPhen() [permulationFunctions.R]
#Modified to set root_sp to first taxa listed in each indivdual gene tree instead of the same root species for every gene tree
#This is to deal with the fact that not every species is in every gene tree due to filtering, missing data, etc
#Elysia says the outgroup/root species is unimportant for permulations, so it is okay to set a different root_sp for each gene tree, and to set the outgroup arbitrarily this way
#Also runs gtSistersList() for SSM permulations to generate a new sisters list for each gene tree (this accounts for missing taxa in each gene tree, which can result in sisters pairings of foreground clades different from that in the master tree)

#'Produces binary permulations for a gene
#' @param tree Tree of the gene of interest (if permmode="cc", set this as the masterTree in trees (i.e., the output from \code{\link{readTrees}}))
#' @param numperms An integer number of permulations
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @return output.list a list containing the set of binary permulated trees
#' @export
generatePermulatedBinPhenQuiet=function(tree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc"){
  if (permmode=="cc"){
    tree_rep = lapply(1:numperms,rep_tree,tree=trees)
    #THIS IS THE PART I ADDED (EK, 11/21/2022) - arbitrarily sets root_sp to the first tip label in the tree
    root_sp=tree_rep[[1]]$tip.label[1]
    permulated.binphens = lapply(tree_rep, simBinPhenoCC,mastertree=trees$masterTree,root_sp=root_sp, fg_vec=fg_vec,sisters_list=sisters_list,pathvec=pathvec,plotTreeBool=F)
  } else if (permmode=="ssm"){
    tree_rep = lapply(1:numperms,rep_tree,tree=tree)
    #THIS IS THE PART I ADDED (EK, 11/21/2022)
    #arbitrarily sets root_sp to the first tip label in the tree
    #root_sp=tree_rep[[1]]$tip.label[1]
    #randomly root the tree
    root_sp=sample(tree_rep[[1]]$tip.label, 1)
    #print(tree_rep[[1]])
    #print(tree_rep[[1]]$tip.label)
    #print(paste("ROOT SPECIES:",root_sp))
    #gets sisters list for each gene tree, to account for missing data and possible different sisters structure in each gene tree
    temp_sis=gtSistersList(tree, fg_vec)
    #OLD SSM version - based on each gene tree
    #permulated.binphens = lapply(tree_rep,simBinPhenoSSM,trees=trees,root_sp=root_sp,fg_vec=fg_vec,sisters_list=temp_sis,pathvec=pathvec, plotTreeBool=T)
    #NEW SSM version - based on master tree but subsets to only include species in gene tree
    permulated.binphens = lapply(tree_rep,simBinPhenoSSM_fromMasterTree,trees=trees,fg_vec=fg_vec,sisters_list=temp_sis,pathvec=pathvec, plotTreeBool=F)
  } else if (permmode=="rank"){
    tree_rep = lapply(1:numperms,rep_tree,tree=trees)
    root_sp=tree_rep[[1]]$tip.label[1]
    permulated.binphens = lapply(tree_rep,simBinPhenoRankQuiet,root_sp=root_sp,fg_vec=fg_vec,sisters_list=sisters_list,plotTreeBool=F)
  } else {
    stop("Invalid binary permulation mode.")
  }
  output.list <- list()
  output.list[[1]] <- permulated.binphens
  return(output.list)
}
#' @keywords internal
rep_tree = function(num_input,tree){
  return(tree)
}

#EK - removed the requirement that permulated tree structure matches true structure for the gene tree (commented out testcondition in bottom while loop)
#EK - added a counter to the while loops such that it only tries 50 times to find a permulated tree the matches the conditions (same number of foreground branches); if it cannot find one after 50 tries it returns a NULL tree
#'Produces one SSM binary permulation for a gene
#' @param tree Tree of the gene of interest
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @return A SSM binary permulated tree
#' @export
simBinPhenoSSM=function(tree, trees, root_sp, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  require(phytools)
  #print("In simBinPhenoSSM")
  #print(tree)
  tip.labels = tree$tip.label # the set of species that exist in the gene tree
  #print(tip.labels)
  ind_fg = which(tip.labels %in% fg_vec) # indices of the observed foreground animals that exist in the gene tree
  #print(ind_fg)

  if (length(ind_fg) == 0){
   #if (length(ind_fg) < 3){ #Seems to stall out later if only 2 foregrounds; trying this to see if it avoids that problem
    t = tree
    t$edge = NULL
    t$edge.length = NULL
    t$Nnode = NULL
    t$tip.label = NULL
  } else {
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree
    #print(paste("Length of foreground tips for this tree:", length(fg_k)))
    #print(fg_k)
    #print(sisters_list)
    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
    #print("Got foreground clades")
    fg_tree = res$tree
    fg.table = res$fg.sisters.table

    print("ROOT SPECIES:")
    print(root_sp)
    #Chronos should make the tree ultrametric to avoid long branch bias during simulations
    #sim.t=chronos(root.phylo(tree, root_sp, resolve.root = T))
    #force.ultrametric forces the tree to be ultrametric but without proper statistical methods
    #sim.t=root.phylo(force.ultrametric(tree, method="nnls"), root_sp, resolve.root=T)
    #Change all branch lengths to 1
    #sim.t=compute.brlen(root.phylo(tree, root_sp, resolve.root = T), 1)
    #Make tree ultrametric AFTER setting branch lengths to 1
    #sim.t=root.phylo(force.ultrametric(t, method="extend"), root_sp, resolve.root=T)
    #Midpoint root the tree
    #sim.t=midpoint.root(tree)
    #Midpoint root AND make tree ultrametric by extending
    #sim.t=force.ultrametric(midpoint.root(tree), method="extend")
    #Midpoint root, branch lengths to 1, make ultrametric
    #sim.t=force.ultrametric(compute.brlen(midpoint.root(tree), 1), method="extend")
    #Original - just roots the tree w/o any ultrametric or branch length modifications
    t=root.phylo(tree, root_sp, resolve.root = T)
    #print(sim.t)
    #print(sim.t$edge.length)
    #print(is.binary(sim.t))
    #write.tree(sim.t, "midpoint_root.bl1.forceUlt_extend.10loci.tre", append=TRUE)
    rm=ratematrix(t, pathvec)
    #Univariate q matrix for testing - MAY NEED TO CHANGE THIS TO REFLECT REAL DATA
    #rm=list(rbind(c(-.75, .75), c(.75, -.75)))
    #Hard code a coonstant rate matrix variance value
    #rm=0.02
    #print(rm)

    if (!is.null(sisters_list)){
      fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
      #print("Got binary perm inputs")
      num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
      num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
      num_tip_sisters_true = length(num_tip_sisters_true)
      #fg_tree_depth_order = getDepthOrder(fg_tree)
    }
    
    #I think this needs to be out of the if(!is.null) because later fg_tree_depth_order needs to exist regardless of if sisters_list is null of not
    fg_tree_depth_order = getDepthOrder(fg_tree)
    #print("Got depth order")
    fgnum = length(which(fg_tree$edge.length == 1))
    if (!is.null(sisters_list)){
      internal = nrow(fg.table)
    } else {
      internal = 0
    }
    tips=fgnum-internal # the number of tips
    #print(paste("Number of foregrounds:", fgnum))
    testcondition=FALSE
    while(!testcondition){ 
      blsum=0
      #Modified so this tries 50 times to generate a tree with the same number of foregrounds; if it doesn't in that time, move on
      #Based on a previous run, <1% of gene trees took 50 times or more, so this shouldn't get rid of that many
      try_count=0
      #while(blsum!=fgnum){
      while( (blsum!=fgnum) && (try_count < 50) ){
        sims=sim.char(t, rm, nsim = 1, model="BM")
        #sims=sim.char(sim.t, rm, nsim = 1, model="BM")
        #sims=sim.char(t, rm, model="discrete", nsim = 1)
        #Simulate evolution of discrete character across phylogeny
        #sims=rTraitDisc(t, model="SYM")
	print(sims)
        nam=rownames(sims)
        s=as.data.frame(sims)
        simulatedvec=s[,1]
        names(simulatedvec)=nam
        top.all=names(sort(simulatedvec, decreasing = TRUE))
        #top.all=names(simulatedvec)[which(simulatedvec==2)]
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
	#print(tree$tip.label)
        #print(length(fg_k))
        #top = sample(tree$tip.label, length(fg_k))
        print("Top:")
	print(top)
        t=foreground2Tree(top, trees, clade="all", plotTree = F, useSpecies=tip.labels)
        blsum=sum(t$edge.length)
	try_count=try_count+1
	print(paste("blsum and try_count:", blsum, try_count))
      }
      if(try_count==50){
        #print("Assigning null tree")
        t = tree
        t$edge = NULL
        t$edge.length = NULL
        t$Nnode = NULL
        t$tip.label = NULL
        testcondition=TRUE
        print("NULL TREE")
      } else{
        t_info = getBinaryPermulationInputsFromTree(t)
        if (!is.null(sisters_list)){
    	  #print("In testconidtion")
          num_tip_sisters_fake = unlist(t_info$sisters_list)
          num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
          num_tip_sisters_fake = length(num_tip_sisters_fake)
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
          #  (num_tip_sisters_fake == num_tip_sisters_true)
	  testcondition = TRUE
        } else {
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
	  testcondition = TRUE
        }
      }  
    }
  }
  #print("Done with simBinPhenoSSM")
  if (plotTreeBool){
    print("HERE")
    if(!(is.null(t$tip.label))){
        print(t)
        plot(t)
        write.tree(t, "temp.tre", append=T)
    } else{
        write("NULL", "temp.tre", append=T)
    }
  }
  return(t)
}

#EK - fixed variable naming bug; "t" is still the output of foreground2Tree and "sim.t" is now the midpoint rooted tree for simulating trait values across
#EK - hard-coded variance as 0.02 for simulations, instead of trying to calculate variance for a simulated continuous trait based on binary data
#EK - subset taxa from the master tree based on the gene tree, instead of giving it the real gene tree, to account for the fact that many branches have length zero in the gene trees
#EK - removed the requirement that permulated tree structure matches true structure for the gene tree (commented out testcondition in bottom while loop)
#EK - added a counter to the while loops such that it only tries 50 times to find a permulated tree the matches the conditions (same number of foreground branches); if it cannot find one after 50 tries it returns a NULL tree
#'Produces one SSM binary permulation for a gene
#' @param tree Tree of the gene of interest
#' @param trees treesObj from \code{\link{readTrees}}
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @return A SSM binary permulated tree
#' @export
simBinPhenoSSM_fromMasterTree=function(tree, trees, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  require(phytools)
  tip.labels = tree$tip.label # the set of species that exist in the gene tree
  ind_fg = which(tip.labels %in% fg_vec) # indices of the observed foreground animals that exist in the gene tree

  if (length(ind_fg) == 0){
    t = tree
    t$edge = NULL
    t$edge.length = NULL
    t$Nnode = NULL
    t$tip.label = NULL
  } else {
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree
    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
    fg_tree = res$tree
    fg.table = res$fg.sisters.table

    #EK - keeps the topology and branch lengths of the master tree but ONLY keeps the tips present in this gene tree
    #EK - also midpoint roots the tree
    #EK - variable is named "sim.t" to distinguish it from "t" (output of foreground2Trees below)
    sim.t=midpoint.root(keep.tip(trees$masterTree, tip.labels))
    #print(sim.t)
    #write.tree(sim.t, "midpoint_root.subset_masterTree.10loci.tre", append=TRUE)
    #EK - hard code a constant rate matrix variance value
    rm=0.02
    #print(rm)

    if (!is.null(sisters_list)){
      fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
      num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
      num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
      num_tip_sisters_true = length(num_tip_sisters_true)
      #fg_tree_depth_order = getDepthOrder(fg_tree)
    }

    #EK - I think this needs to be out of the if(!is.null) because later fg_tree_depth_order needs to exist regardless of if sisters_list is null of not
    fg_tree_depth_order = getDepthOrder(fg_tree)
    fgnum = length(which(fg_tree$edge.length == 1))
    if (!is.null(sisters_list)){
      internal = nrow(fg.table)
    } else {
      internal = 0
    }
    tips=fgnum-internal # the number of tips
    #print(paste("Number of foregrounds:", fgnum))
    testcondition=FALSE
    while(!testcondition){
      blsum=0
      #EK - Modified so this tries 50 times to generate a tree with the same number of foregrounds; if it doesn't in that time, move on
      #EK - Based on a previous run, <1% of gene trees took 50 times or more, so this shouldn't get rid of that many
      try_count=0
      #while(blsum!=fgnum){
      while( (blsum!=fgnum) && (try_count < 50) ){
        sims=sim.char(sim.t, rm, nsim = 1, model="BM")
        #print(sims)
        nam=rownames(sims)
        s=as.data.frame(sims)
        simulatedvec=s[,1]
        names(simulatedvec)=nam
        top.all=names(sort(simulatedvec, decreasing = TRUE))
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
        #print("Top:")
        #print(top)
        t=foreground2Tree(top, trees, clade="all", plotTree = F, useSpecies=tip.labels)
        blsum=sum(t$edge.length)
        try_count=try_count+1
        #print(paste("blsum and try_count:", blsum, try_count))
      }
      if(try_count==50){
        #print("Assigning null tree")
        t = tree
        t$edge = NULL
        t$edge.length = NULL
        t$Nnode = NULL
        t$tip.label = NULL
        testcondition=TRUE
        #print("NULL TREE")
      } else{
        t_info = getBinaryPermulationInputsFromTree(t)
        if (!is.null(sisters_list)){
          num_tip_sisters_fake = unlist(t_info$sisters_list)
num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
          num_tip_sisters_fake = length(num_tip_sisters_fake)
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
          #  (num_tip_sisters_fake == num_tip_sisters_true)
          testcondition = TRUE
        } else {
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
          testcondition = TRUE
        }
      }
    }
  }
  if (plotTreeBool){
    if(!(is.null(t$tip.label))){
        print(t)
        plot(t)
        write.tree(t, "temp.tre", append=T)
    } else{
        write("NULL", "temp.tre", append=T)
    }
  }
  return(t)
}

####generatePermulatedBinPhenSSMBatchedQuiet()####
#EK - calls my modified version of generatePermulatedBinPhen

#'Produces binary SSM permulations for a list of genes
#' @param trees_list A list containing the trees of all genes of interest (formatted like trees in treesObj from \code{\link{readTrees}})
#' @param numperms An integer number of permulations
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @return simPhenoList A list containing binary permulated trees for each gene
#' @export
generatePermulatedBinPhenSSMBatchedQuiet=function(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec){
  masterTree = trees$masterTree
  master.tips = masterTree$tip.label
  df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=masterTree,foregrounds=fg_vec)
  df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
  attr = attributes(df.list[[1]])
  col_names = attr$names
  attr2 = attributes(df.list)
  row_names = attr2$names

  colnames(df.converted) = col_names
  rownames(df.converted) = row_names

  df.converted$num.fg = as.integer(df.converted$num.fg)
  df.converted$num.spec = as.integer(df.converted$num.spec)

  spec.members = df.converted$spec.members

  # Group gene trees based on the similarity of their species membership
  grouped.trees = groupTrees(spec.members)
  ind.unique.trees = grouped.trees$ind.unique.trees
  ind.unique.trees = unlist(ind.unique.trees)
  ind.tree.groups = grouped.trees$ind.tree.groups

  # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
  unique.trees = trees_list[ind.unique.trees]

  # Generate simulated phenotypes
  unique.pheno.list = mapply(generatePermulatedBinPhenQuiet,unique.trees,MoreArgs = list(numperms=numperms,trees=trees,root_sp=root_sp,fg_vec=fg_vec,sisters_list=sisters_list,pathvec=pathvec,permmode="ssm"))
  # Allocate the simulated phenotypes for unique trees to their respective groups
  simPhenoList = vector("list", length = length(trees_list))
  for (j in 1:length(simPhenoList)){
    simPhenoList[[j]] = vector("list",length=numperms)
  }
  for (i in 1:length(unique.pheno.list)){
    ind.unique.tree = ind.unique.trees[i]
    ind.tree.group = ind.tree.groups[[i]]
    unique.pheno = unique.pheno.list[[i]]
    for (k in 1:length(ind.tree.group)){
      simPhenoList[[ind.tree.group[k]]] = unique.pheno
    }
  }

  attributes(simPhenoList)$names = row_names

  return(simPhenoList)
}


#NOTE: Below are some functions I didn't change but copied because they are marked as "internal" so the functions that call them can't find them unless they're in the same file

#' @keywords internal
getSpeciesMembershipStats = function(tree,masterTree,foregrounds){
  master.tips = masterTree$tip.label
  tips = tree$tip.label
  spec_membership = which(master.tips %in% tips)
  fg_membership = which(foregrounds %in% tips)

  num_spec = length(spec_membership)
  num_fg = length(fg_membership)

  spec.members = rep(0,length(master.tips))
  spec.members[spec_membership] = 1
  spec.members = toString(spec.members)

  fg.members = rep(0,length(foregrounds))
  fg.members[fg_membership] = 1
  fg.members = toString(fg.members)

  df = data.frame("num.fg"=as.character(num_fg), "num.spec"=as.character(num_spec), "spec.members"=spec.members, "fg.members"=fg.members)

  return(df)
}

#' @keywords internal
groupTrees = function(spec.members){
  unique.trees = unique(spec.members)
  ind.tree.groups = lapply(unique.trees,findGroupedTrees,spec.members=spec.members)
  ind.unique.trees = lapply(ind.tree.groups, function(i) min(i))
  output.list = list()
  output.list$ind.unique.trees = ind.unique.trees
  output.list$ind.tree.groups = ind.tree.groups
  return(output.list)
}

#' @keywords internal
findGroupedTrees = function(unique.tree,spec.members){
  ind.grouped.trees = which(spec.members == unique.tree)
  return(ind.grouped.trees)
}

#'Calculates the clade mappings between the gene tree and the master tree (with the complete topology)
#' @param gene_tree A binary phenotype tree of a gene
#' @param treesObj treesObj from \code{\link{readTrees}}
#' @return output.map A list containing a dataframe of clades mapping
#' @export
matchAllNodesClades=function(gene_tree, treesObj){
  foregrounds = getForegroundsFromBinaryTree(gene_tree)
  tree1 = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F)

  map=matchNodesInject_mod(tree1,treesObj$masterTree)
  map=map[order(map[,1]),]
  #map

  output.map = list()
  output.map[[1]]=map
  output.map
}

#' @keywords internal
compareClades=function(clade2index,desc.tree2,clade1){
  output=NA
  clade2 = desc.tree2[[clade2index]]
  if (all(clade1 %in% clade2)){
    output = as.numeric(names(desc.tree2)[clade2index])
  }
  return(output)
}

#' @keywords internal
findMappedClade2Node=function(desc.tr1.index,desc.tr2.index.list,desc.tree1,desc.tree2){
  clade1 = desc.tree1[[desc.tr1.index]]
  mapped.clades.list = lapply(desc.tr2.index.list,compareClades,desc.tree2=desc.tree2,clade1=clade1)
  mapped.clades = unlist(mapped.clades.list)
  mapped.clades.nonNA = mapped.clades[!is.na(mapped.clades)]
  return(max(mapped.clades.nonNA))
}

#' @keywords  internal
matchNodesInject_mod=function (tr1, tr2){
  if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
    #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
    stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
  }
  commontiplabels <- intersect(tr1$tip,tr2$tip)
  if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
    stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
  }
  #if(RF.dist(tr1,tr2)>0){
  #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
  #}

  toRm=setdiff(tr2$tip.label, tr1$tip.label)
  desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                              x)$tip.label)
  names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
  desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                              x)$tip.label)
  names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
  Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                           c("tr1", "tr2")))
  Nodes[,1] = as.numeric(names(desc.tr1))
  desc.tr1.index.list = as.list(1:length(desc.tr1))
  desc.tr2.index.list = as.list(1:length(desc.tr2))


  mapped.clade.list = lapply(desc.tr1.index.list,findMappedClade2Node,desc.tr2.index.list,desc.tree1=desc.tr1,desc.tree2=desc.tr2)
  Nodes[,2] = unlist(mapped.clade.list)

  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  if(any(table(Nodes[,2])>1)){
    stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
  }

  Nodes
}

#' @keywords internal
calculatePermulatedPaths=function(permulated.trees,map,treesObj){
  permulated.paths=lapply(permulated.trees,tree2PathsCladesQuiet,trees=treesObj)
  output.list = list()
  output.list[[1]] = permulated.paths
  output.list
}

#EK - Changed to call quiet version of tree2Paths_map
#'A modification of the tree2Paths function that takes in pre-calculated mappings
#' @param tree the input tree to be converted into paths
#' @param trees treesObj from \code{\link{readTrees}}
#' @export
tree2PathsCladesQuiet=function(tree,trees){
  map = matchAllNodesClades(tree,trees)
  path = tree2Paths_mapQuiet(tree,map[[1]],trees)
  names(path) = colnames(trees$paths)
  path
}

#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' @keywords internal
tree2Paths_mapQuiet=function(tree, map, treesObj, binarize=NULL, useSpecies=NULL){
  if (class(tree)[1]=="phylo"){
    stopifnot(class(tree)[1]=="phylo")
    stopifnot(class(treesObj)[2]=="treesObj")

    if (is.null(tree$tip.label)){
      vals=as.double(rep(NA,length(treesObj$ap$dist)))
    } else {
      foregrounds = getForegroundsFromBinaryTree(tree)
      tree = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F)


      isbinarypheno <- sum(tree$edge.length %in% c(0,1)) == length(tree$edge.length) #Is the phenotype tree binary or continuous?
      if (is.null(binarize)) { #unless specified, determine default for binarize based on type of phenotype tree
        if (isbinarypheno) {
          binarize = T #default for binary phenotype trees: set all positive paths = 1
        } else {
          binarize = F #default for continuous phenotype trees: do not convert to binary
        }
      }

      #unroot if rooted
      if (is.rooted(tree)) {
        tree = unroot(tree)
      }

      #reduce tree to species in master tree and useSpecies
      sp.miss = setdiff(tree$tip.label, union(treesObj$masterTree$tip.label, useSpecies))
      #if (length(sp.miss) > 0) {
      #  message(paste0("Species from tree not present in master tree or useSpecies: ", paste(sp.miss,
      #                                                                                       collapse = ",")))
      #}

      if (!is.null(useSpecies)) {
        tree = pruneTree(tree, intersect(intersect(tree$tip.label, treesObj$masterTree$tip.label), useSpecies))
      } else {
        tree = pruneTree(tree, intersect(tree$tip.label, treesObj$masterTree$tip.label))
      }
      treePaths=allPaths(tree)

      #remap the nodes
      treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
      treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]

      #indices for which paths to return
      ii=treesObj$ap$matIndex[(treePaths$nodeId[,2]-1)*nrow(treesObj$ap$matIndex)+treePaths$nodeId[,1]]

      vals=double(length(treesObj$ap$dist))
      vals[]=NA
      vals[ii]=treePaths$dist
      if(binarize){
        if(isbinarypheno) {
          vals[vals>0]=1
        } else {
          mm=mean(vals)
          vals[vals>mm]=1
          vals[vals<=mm]=0
        }
      }
    }
  } else {
    vals=as.double(rep(NA,length(treesObj$ap$dist)))
  }

  vals
}

#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' Generate a phenotype paths vector from a phenotype tree
#'
#' \code{tree2Paths} generates a phenotype paths vector matching the treesObject
#'     from a tree where branches specify phenotypes.
#'
#' The tree topology of the phenotype tree must match that of the master tree within the treesObject.
#'
#' @param tree A phenotype tree, with branch length encoding a phenotype.
#' @param treesObj A treesObject created by \code{\link{readTrees}}
#' @param binarize Force binary path representation. Default action depends upon the type of data within the phenotype tree
#'     (binary or continuous).
#'     \itemize{
#'    \item If binary (all branch lengths == 0 or 1): Sets all positive path values to 1. Useful if the tree has non-zero branch lengths
#'        for an internal branch or branches; otherwise, values are simply added along branches when calculating paths.
#'        Default behavior: binarize = TRUE.
#'    \item If continuous (not all branch lengths == 0 or 1): Sets all path values > the mean to 1 and all those <= the mean to 0.
#'        Converts a continuous phenotype to a binary phenotype, with state determined by comparison to the mean across all paths.
#'        Default behavior: binarize = FALSE.
#'        }
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A vector of length equal to the number of paths in treesObj
#' @export
tree2Paths=function(tree, treesObj, binarize=NULL, useSpecies=NULL, categorical = F){
  stopifnot(class(tree)[1]=="phylo")
  stopifnot(class(treesObj)[2]=="treesObj")

  isbinarypheno <- sum(tree$edge.length %in% c(0,1)) == length(tree$edge.length) #Is the phenotype tree binary or continuous?
  if (is.null(binarize)) { #unless specified, determine default for binarize based on type of phenotype tree
    if (isbinarypheno) {
      binarize = T #default for binary phenotype trees: set all positive paths = 1
    } else {
      binarize = F #default for continuous phenotype trees: do not convert to binary
    }
  }
  #unroot if rooted
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }

  #reduce tree to species in master tree and useSpecies
  sp.miss = setdiff(tree$tip.label, union(treesObj$masterTree$tip.label, useSpecies))
  if (length(sp.miss) > 0) {
    #message(paste0("Species from tree not present in master tree or useSpecies: ", paste(sp.miss,
    #                                                                                     collapse = ",")))
  }
  if (!is.null(useSpecies)) {
    tree = pruneTree(tree, intersect(intersect(tree$tip.label, treesObj$masterTree$tip.label), useSpecies))
  } else {
    tree = pruneTree(tree, intersect(tree$tip.label, treesObj$masterTree$tip.label))
  }

  treePaths=allPaths(tree, categorical = categorical)
  map=matchAllNodes(tree,treesObj$masterTree)

  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]

  #indices for which paths to return
  ii=treesObj$ap$matIndex[(treePaths$nodeId[,2]-1)*nrow(treesObj$ap$matIndex)+treePaths$nodeId[,1]]

  vals=double(length(treesObj$ap$dist))
  vals[]=NA
  vals[ii]=treePaths$dist
  if(binarize){
    if(isbinarypheno) {
      vals[vals>0]=1
    } else {
      mm=mean(vals)
      vals[vals>mm]=1
      vals[vals<=mm]=0
    }
  }
  vals
}



#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' Creates a binary trait tree from a set of foreground species.
#' @param foreground. A character vector containing the foreground species
#' @param treesObj A treesObj created by \code{\link{readTrees}}
#' @param collapse2anc Put all the weight on the ancestral branch when the trait appears on a whole clade
#' (redundant to "clade", kept for backwards compatibility)
#' @param plotTree Plot a tree representation of the result
#' @param wholeClade Whether to implement the weighted edge option across
#' all members of a foreground clade (redundant to "clade", kept for backwards compatibility)
#' @param clade A character string indicating which branches within the clade
#' containing the foreground species should be set to foreground. Must be one
#' of the strings "ancestral", "terminal", "all".
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @param weighted if set to TRUE weights foreground edges belonging to the same clade such that their branch lengths sum up to 1 (only done for clade options "all" and "terminal").
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A tree with edge.lengths representing phenotypic states
#' @export
foreground2Tree = function(foreground,treesObj, plotTree=T, clade=c("ancestral","terminal","all"), weighted = F, transition = "unidirectional", useSpecies=NULL){
  clade <- match.arg(clade)
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      #message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
      #                                                                            collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  } else {
    useSpecies = res$tip.label
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0,length(res$edge.length))
  if(clade == "terminal"){
    res$edge.length[nameEdges(res) %in% foreground] = 1
    names(res$edge.length) = nameEdges(res)
  }else if(clade == 'ancestral'){
    weighted = F
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }
  }else{
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }
  }
  if(weighted){
    if(clade == 'all'){
      tobeweighted <- rep(TRUE,length(res$edge.length))
      tobeweighted[res$edge.length == 0] <- FALSE
      while(sum(tobeweighted)>0){
        edgetodo <- which(tobeweighted == T)[1]
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(length(clade.down.edges) > 1){
          clade.edges = c(clade.down.edges, edgetodo)
          clade.edges.toweight <- clade.edges[res$edge.length[clade.edges] == 1]
          res$edge.length[clade.edges.toweight] <- 1.0/(length(clade.edges.toweight))
          tobeweighted[clade.edges] <- FALSE
        } else{
          tobeweighted[clade.down.edges] <- FALSE
        }
      }
    } else if(clade == 'terminal'){
      tobeweightededgeterminalnode <- unique(res$edge[(res$edge[,2] %in% c(1:length(res$tip.label))),1])
      tobeweighted <- setdiff(match(tobeweightededgeterminalnode, res$edge[,2]), NA)
      for(edgetodo in tobeweighted){
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(all(res$edge.length[clade.down.edges]==1)){
          res$edge.length[clade.down.edges] <- 0.5
        }
      }
    }
  }
  if(plotTree){
    res2=res
    mm=min(res2$edge.length[res2$edge.length>0])
    res2$edge.length[res2$edge.length==0]=max(0.02,mm/20)
    plot(res2, main = paste0("Clade: ",clade,'\nTransition: ',transition,'\nWeighted: ',weighted), cex = 0.5)
    if(weighted){
      labs <- round(res$edge.length,3)
      labs[labs == 0] <- NA
      edgelabels(labs, col = 'black', bg = 'transparent', adj = c(0.5,-0.5),cex = 0.4,frame='n')
    }
  }
  res
}

#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' main RER computation function
#' @param treesObj A treesObj created by \code{\link{readTrees}}
#' @param a cutoff value for branch lengths bellow which the branch lengths will be discarded, very data dependent but should roughly correspond to 0 or 1 sequence change on that branch. If left NULL this whill be set to the bottom 0.05 quantile. Set to 0 for no cutoff.
#' @param transform The transformation to apply to the trees branch values before computing relative rates. Available options are sqrt and log, sqrt is recommended.
#' @param weighted Use weighted regression to compute relative rates, meant to correct for the non-constant mean-variance relationship in evolutionary rate data.
#' @param useSpecies Give only a subset of the species to use for RER calculation. Some times excluding unusually long branches can provide more stable results
#' @param min.sp The minimum number of species needed to compute RER
#' @param scale Scale relative rates internally for each species subset. Increases computation time with little apparent benefit. Better to scale the final matrix.
#' @param doOnly The index of a specific tree in the treesObj to calculate RER for. Useful if a single result is needed quickly.
#' @param maxT The maximum number of trees to compute results for. Since this function takes some time this is useful for debugging.
#' @param plot Whether to plot the output of the correction for mean-variance relationship.
#' @return A numer of trees by number of paths matrix of relative evolutionary rates. Only an independent set of paths has non-NA values for each tree.
#' @export
getAllResiduals=function(treesObj, cutoff=NULL, transform="sqrt", weighted=T,  useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05, plot=T){

  if(is.null(cutoff)){
    cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
    message(paste("cutoff is set to", cutoff))
  }
  if (weighted){
    weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=plot)
    residfunc=fastLmResidMatWeighted
  }
  else{
    residfunc=fastLmResidMat
  }
  # residfunc=naresid

  if (is.null(useSpecies)){
    useSpecies=treesObj$masterTree$tip.label
    #mappedEdges=trees$mappedEdges
  }
  if(is.null(maxT)){
    maxT=treesObj$numTrees
  }
  if(transform!="none"){
    transform=match.arg(transform,c("sqrt", "log"))
    transform=get(transform)
  }
  else{
    transform=NULL
  }



  #cm is the names of species that are included in useSpecies and the master tree
  cm=intersect(treesObj$masterTree$tip.label, useSpecies)
  sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
  if (length(sp.miss) > 0) {
    #message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
    #                                                                             collapse = ",")))

  }

  rr=matrix(nrow=nrow(treesObj$paths), ncol=ncol(treesObj$paths))

  #maximum number of present species
  maxn=rowSums(treesObj$report[,cm])

  if(is.null(doOnly)){
    doOnly=1
  }
  else{
    maxT=1
  }
  skipped=double(nrow(rr))
  skipped[]=0

  for (i in doOnly:(doOnly+maxT-1)){

    if(sum(!is.na(rr[i,]))==0&&!skipped[i]==1){


      #get the ith tree
      tree1=treesObj$trees[[i]]

      #get the common species, prune and unroot
      both=intersect(tree1$tip.label, cm)
      if(length(both)<min.sp){
        next
      }
      tree1=unroot(pruneTree(tree1,both))

      #do the same for the refTree


      #find all the genes that contain all of the species in tree1
      allreport=treesObj$report[,both]
      ss=rowSums(allreport)
      iiboth=which(ss==length(both)) #this needs to be >1
      if (length(iiboth) < 2) {
        message(paste("Skipping i =",i,"(no other genes with same species set)"))
        next
      }

      nb=length(both)
      ai=which(maxn[iiboth]==nb)


      message(paste("i=", i))


      if(T){

        ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)

        ii= treesObj$matIndex[ee[, c(2,1)]]

        allbranch=treesObj$paths[iiboth,ii]
        if (is.null(dim(allbranch))) {
          message(paste("Issue with gettiing paths for genes with same species as tree",i))
          return(list("iiboth"=iiboth,"ii"=ii))
        }

        if(weighted){
          allbranchw=weights[iiboth,ii]
        }
        if(scaleForPproj){
          nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
        }
        else{
          nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
        }

        iibad=which(allbranch<cutoff)
        #don't scale
        #allbranch=scaleMat(allbranch)
        if(!is.null(transform)){
          nv=transform(nv)
          allbranch=transform(allbranch)
        }
        allbranch[iibad]=NA




        if(!scale){
          if(!weighted){
            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv))

          }
          else{

            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv), allbranchw[ai, ,drop=F])

          }
        }

        else{

          if(!weighted){
            proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv))
          }
          else{

            proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv),allbranchw)
          }

          proj=scale(proj, center = F)[ai, , drop=F]

        }


        #we have the projection



        rr[iiboth[ai],ii]=proj

      }

    }}
  message("Naming rows and columns of RER matrix")
  rownames(rr)=names(treesObj$trees)
  colnames(rr)=namePathsWSpecies(treesObj$masterTree)
  rr
}


#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' Creates a categorical trait tree from a set of tip species.
#'@param tipvals the trait/phenotype/character value at the tip, \code{names(tip.vals)} should match some of the \code{mastertree$tip.label}, though a perfect match is not required
#'@param treesObj A treesObj created by \code{\link{readTrees}}
#'@param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#'@param model Specifies what rate model to use
#'@param plot Plots a phenotype tree
#'@param anctrait The trait to use for all ancestral species instead of inferring ancestral states if not NULL. The default is NULL.
#'@return A tree with edge.length representing phenotype states
#'@export
char2TreeCategorical <- function (tipvals, treesObj, useSpecies = NULL,
                                     model = "ER", root_prior = "auto",
                                     plot = FALSE, anctrait = NULL)
{
  # get the master tree and prune to include useSpecies/species with phenotype data
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      #message(paste0("Species from master tree not present in useSpecies: ",
      #               paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label,
                                                 names(tipvals)))
  }
  # use ASR to infer phenotype tree
  if (is.null(anctrait)) {

    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)

    ancliks = getAncLiks(mastertree, intlabels$mapped_states, rate_model = model,
                         root_prior = root_prior)

    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i,])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]

    # convert to binary tree if necessary, plot, & return tree
    if(length(unique(tipvals)) == 2) {
      if(sum(! unique(tipvals) %in% c(TRUE,FALSE)) > 0) { # check that the two categories are TRUE/FALSE
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      } else {
        tree$edge.length = ifelse(tree$edge.length == 2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical(tree, category_names = intlabels$state_names,
                          master = mastertree, node_states = states)
    }
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot,
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if(plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals ==
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        # get states for plotting
        states = res$edge.length[order(res$edge[,2])]
        states = c(j, states) # add root since it's not included in res$edge.length (no edge leading to the root)
        plotTreeCategorical(res, category_names = traits,
                            master = treesObj$masterTree,
                            node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}

#' @keywords internal
inferUnidirectionalForegroundClades <- function(tree, fgd = NULL, ancestralOnly = F){
  finaltree <- tree
  finaltree$edge.length <- rep(0, length(tree$edge.length))
    finaltree$edge.length[nameEdges(finaltree) %in% fgd] <- 1
  #figure out node depth - terminal nodes have depth of 1; higher numbers indicate ancestral nodes;
  nodedepths <- node.depth(finaltree)
  edgeterminalnodedepths <- nodedepths[finaltree$edge[,2]]
  #going from 1-away from terminal ancestral branch to the base of the tree, figure out branches where all downstream clades are foreground
  for(inode in sort(unique(edgeterminalnodedepths))[-1]){
    edgesToDo <- which(edgeterminalnodedepths == inode)
    for(edgeindex in edgesToDo){
      clade.edges = getAllCladeEdges(finaltree, edgeindex)
      if(all(finaltree$edge.length[clade.edges]==1)){
        finaltree$edge.length[edgeindex] <- 1
      }
    }
  }
  if(ancestralOnly){
    for(edgeii in 1:length(finaltree$edge.length)){
      if(finaltree$edge.length[edgeii] == 1){
        if(nameEdges(finaltree)[edgeii]==""){
          clade.edges = setdiff(getAllCladeEdges(finaltree, edgeii), edgeii)
          finaltree$edge.length[clade.edges] <- 0
        }
      }
    }
  }
  finaltree
}

#' @keywords internal
nameEdges=function(tree){
  nn=character(nrow(tree$edge))
  iim=match(1:length(tree$tip.label), tree$edge[,2])
  nn[iim]=tree$tip.label
  nn
}

getAllCladeEdges=function(tree, AncEdge){
  node=tree$edge[AncEdge,2]
  #get descendants
  iid=getDescendants(tree, node)
  #find their edges
  iim=match(iid, tree$edge[,2])
  iim
}

#' @keywords internal
getDepthOrder=function(fgTree){
  unq_edge_lengths = unique(fgTree$edge.length)
  if (length(which(!(unq_edge_lengths %in% c(0,1)))) > 0){
    stop('Phenotype must be binary.')
  }
  all_edges = fgTree$edge
  num_tip_species = length(fgTree$tip.label)

  idx_fg_branches = which(fgTree$edge.length == 1)
  if (length(idx_fg_branches)==1){
	fg_edges = fgTree$edge[idx_fg_branches,]
	fg_edges = t(as.data.frame(fg_edges))
  } else {
	fg_edges = fgTree$edge[idx_fg_branches,]
	tip_fg_edges = fg_edges[which(fg_edges[,2] <= num_tip_species),]
	tip_foregrounds = fgTree$tip.label[tip_fg_edges[,2]]
	node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
  }

  idx_node_edges = which(fg_edges[,2] > num_tip_species)
  if (length(idx_node_edges) == 1){
    node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    node_fg_edges = t(as.data.frame(node_fg_edges))
  }
  if (length(idx_node_edges) == 0) {
    sisters_list = NULL
    depth_order = NULL
  } else {
    #node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    daughters_info_list = list()
    parents = NULL
    for (i in 1:nrow(node_fg_edges)){
      edge_i = node_fg_edges[i,]
      # find the daughters of this node
      idx_daughters_i = which(all_edges[,1] == edge_i[2])
      daughter_edges = all_edges[idx_daughters_i,]
      daughters_info_list[[i]] = daughter_edges[,2]
      parents = c(parents, edge_i[2])
    }
    names(daughters_info_list) = parents
    ### write something to order the branches based on depth
    tip_fg_ids = tip_fg_edges[,2]
    depth_order = rep(NA, length(daughters_info_list))
    names(depth_order) = names(daughters_info_list)
    order_assigned = NULL
    while(length(which(is.na(depth_order))) > 0){
      idx_na = which(is.na(depth_order))
      if (length(idx_na) > 0){
        for (j in 1:length(idx_na)){
          idx_na_j = idx_na[j]
          parent_j = parents[idx_na_j]
          daughters_j = daughters_info_list[[idx_na_j]]
          num_tip_daughters = length(which(daughters_j %in% tip_fg_ids))
          if (num_tip_daughters == 2){
            depth_order[idx_na_j] = 1
            order_assigned = c(order_assigned, parent_j)
          } else if (num_tip_daughters==1){
            node_daughter = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (node_daughter %in% order_assigned){
              depth_order[idx_na_j] = depth_order[as.character(node_daughter)] + 1
              order_assigned = c(order_assigned, parent_j)
            }
          } else if (num_tip_daughters==0){
            node_daughters = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (length(which(node_daughters %in% order_assigned)) == 2){
              node_daughters_depths = depth_order[as.character(node_daughters)]
              depth_order[idx_na_j] = max(node_daughters_depths) + 1
              order_assigned = c(order_assigned, parent_j)
            }
          }
        }
      }
    }
  }
  depth_order
}

#' @keywords internal
extractCorResults=function(corResultsList,numperms,mode="Rho"){
  output = matrix(NA,nrow=length(corResultsList),ncol=length(corResultsList[[1]]))
  for (i in 1:length(corResultsList)){
    gene = corResultsList[[i]]
    output[i,] = extractPermulationResults(gene,numperms,mode)
  }
  return(output)
}

#' @keywords internal
extractPermulationResults=function(gene,numperms,mode="Rho"){
  table_perm = lapply(gene,linearizeCorResults)
  df_perm = do.call(rbind,table_perm)
  if (mode=="Rho"){
    output=df_perm[,1]
  } else if (mode=="P"){
    output=df_perm[,3]
  }
  return(output)
}

#' @keywords internal
linearizeCorResults=function(cor_result){
  vec.cor = unlist(cor_result)
  return(vec.cor)
}

#' @keywords internal
allPaths=function(tree){
  dd=dist.nodes(tree) # returns a matrix with col_names and row_names denoting the numbers of the tips and the nodes, containing distances between nodes
  allD=double() # this is the 'path' vector that will be outputted in the end
  nn=matrix(nrow=0, ncol=2) # initialize an empty matrix that will store nodeIds of each node and its ancestors corresponding to each position in the path vector
  nA=length(tree$tip.label)+tree$Nnode # nA is the total number of nodes -- tree$Nnode is the number of internal nodes, length(tree$tip.label) is the number of tip nodes
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1
  # Below here is where the vector of 'paths' are generated. Starting from the first to last tip labels, and then the internal nodes from root to last.
  for ( i in 1:nA){ # for each node i
    ia=getAncestors(tree,i) # find the ancestor nodes of node i
    if(length(ia)>0){ # If node i has ancestors at all
      allD=c(allD, dd[i, ia]) # extends allD per node, where each extension is the distances between node i and its ancestors
      nn=rbind(nn,cbind(rep(i, length(ia)), ia)) # extend the xxx-by-2 matrix by nodeId pairs of node i and its ancestors
      for (j in ia){
        matIndex[i,j]=index
        index=index+1
      }
    }
  }
  return(list(dist=allD, nodeId=nn, matIndex=matIndex))
}

#' @keywords internal
getAncestors=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  if(length(im)==0){
    return()
  }
  else{
    anc=tree$edge[im,1]
    return(c(anc, getAncestors(tree, anc)))
  }
}

#' @keywords internal
findPairs=function(binary.tree){
  tip.labels = binary.tree$tip.label
  edge = binary.tree$edge
  edge.length = binary.tree$edge.length
  ind.fg.edge = which(edge.length == 1)
  fg.edges = edge[ind.fg.edge,]

  # Find the pairs
  fg.pairs.table = matrix(nrow=0,ncol=2)
  colnames(fg.pairs.table) = c("species1","species2")

  if (length(as.vector(fg.edges)) > 2){
    nodeId.start = sort(fg.edges[,1])
    fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
    for (nn in 1:(length(nodeId.start)-1)){
      if (nodeId.start[nn] == nodeId.start[nn+1]){
        fg_ca = c(fg_ca,nodeId.start[nn])
        fg.pairs.table = rbind(fg.pairs.table, fg.edges[which(fg.edges[,1]==nodeId.start[nn]),2])
      }
    }
    rownames(fg.pairs.table) = fg_ca
  }
  fg.pairs.table
}

#'Calculate permulation correlation statistics
#' @param permulated.paths A nested list of permulated paths (e.g., output of \code{\link{calculatePermulatedPaths_apply}}
#' @param RERmat An RER matrix calculated using \code{\link{getAllResiduals}}.
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the foreground (non-zero phenotype values)
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorizetrait=some number and winsorizeRER=some number to override
#' @return A nested list containing the correlation statistics for the permulations
#' @export
calculateCorPermuted=function(permulated.paths,RERmat,min.sp=10,min.pos=2,method="k"){
  corMatList = lapply(permulated.paths,getAllCorSSM,RERmat,min.sp=min.sp,min.pos=min.pos,method=method)
  output.list <- list()
  output.list[[1]] <- corMatList
  return(output.list)
}

#' @keywords internal
getAllCorSSM=function(charP, RERmat, method="auto",min.sp=10, min.pos=2, winsorizeRER=NULL, winsorizetrait=NULL, weighted=F){
  if (method=="auto"){
    lu=length(unique(charP))
    if(lu==2){
      method="k"
      message("Setting method to Kendall")
    }
    else if (lu<=5){
      method="s"
      message("Setting method to Spearman")
    }
    else{
      method="p"
      message("Setting method to Pearson")
      if(is.null(winsorizeRER)){
        message("Setting winsorizeRER=3")
        winsorizeRER=3
      }
      if(is.null(winsorizetrait)){
        message("Setting winsorizetrait=3")
        winsorizetrait=3
      }
    }
  }
  win=function(x,w){
    xs=sort(x[!is.na(x)], decreasing = T)
    xmax=xs[w]
    xmin=xs[length(xs)-w+1]

    x[x>xmax]=xmax
    x[x<xmin]=xmin
    x
  }
  dim(RERmat) <- c(1,length(RERmat))
  corout=matrix(nrow=nrow(RERmat), ncol=3)
  rownames(corout)=rownames(RERmat)

  colnames(corout)=c("Rho", "N", "P")

  for( i in 1:nrow(corout)){

    if(((nb<-sum(ii<-(!is.na(charP)&!is.na(RERmat[i,]))))>=min.sp)){
      if (method!="p"&&sum(charP[ii]!=0)<min.pos){
        next
      }

      if(!weighted){

        x=RERmat[i,]

        #winsorize
        indstouse=which(!is.na(x) & !is.na(charP))
        if(!is.null(winsorizeRER)){
          x=win(x[indstouse], winsorizeRER)
        }else{
          x=x[indstouse]
        }
        if(!is.null(winsorizetrait)){
          y=win(charP[indstouse], winsorizetrait)
        }else{
          y=charP[indstouse]
        }

        cres=cor.test(x, y, method=method, exact=F)
        corout[i,1:3]=c(cres$estimate, nb, cres$p.value)
      }
      else{
        charPb=(charP[ii]>0)+1-1

        weights=charP[ii]
        weights[weights==0]=1

        cres=wtd.cor(RERmat[i,ii], charPb, weight = weights, mean1 = F)
        corout[i, 1:3]=c(cres[1], nb, cres[4])
      }
    }
    else{
      #show(i)
      #show(c(nb, charP[ii]))
    }

  }

  corout=as.data.frame(corout)
  corout$p.adj=p.adjust(corout$P, method="BH")
  corout
}

#EK - new function I wrote
#EK - Make a custom sisters_list for a given gene tree to account for missing data that makes new foreground sisters groups
#EK - As of now (2/14/2023) this function does NOT account for anything more complicated than 2 nested clades within foregrounds; I think I could make this work for an infinite number of nested foreground clades using recursion but don't have time to write that now
gtSistersList=function(tree, fg_vec){
  require(phytools)
  sis_fg<-list()
  for(i in fg_vec){
    #print(paste("Working on foreground", i))
    #Make sure this fg taxon is in the tree
    if(length(getSisters(tree, i)) > 0){
      #Get all sister taxa of i, or all descendents of sister node of i
      sis_tips<-tree$tip.label[getDescendants(tree, getSisters(tree, i))]
      #print(sis_tips)
      #Only add if all sisters/sister node descendents are foreground taxa
      if(all(sis_tips %in% fg_vec)){
        #print(sis_tips)
        #If only one foreground sister, make a new sister vector; append to sis_fg list if this sister pair does not already exist in the sis_fg list
        if(length(sis_tips) == 1){
          if(!(list(sort(c(i, sis_tips))) %in% sis_fg)){
            sis_fg[[paste0("clade",length(sis_fg)+1)]] = sort(c(i, sis_tips))
          }
        } else if(length(sis_tips) == 2){
            #If sis_tips==2, it means there is a sister node with 2 descendents that are also foreground tip; add these sisters as their own clade first, if they are not already in sis_fg, then add this fg taxon plus the clade with the descendents
            if(!(list(sort(sis_tips)) %in% sis_fg)){
              sis_fg[[paste0("clade",length(sis_fg)+1)]] = sort(sis_tips)
            }
            des_index<-match(list(sort(sis_tips)), sis_fg)
            sis_fg[[paste0("clade",length(sis_fg)+1)]] = c(paste0("clade",des_index), i)
        } #else do nothing; if no sisters/descendents are foregrounds we don't want to add anything
          #IMPORTANT NOTE: doesn't account for anything more complicated than 2 nested clades w/ foregrounds
      }
    }
  }
  #print(sis_fg)
  if(length(sis_fg) == 0){
    NULL #Returns NULL if no foreground taxa are sister
  } else{
    sis_fg
  }
}
