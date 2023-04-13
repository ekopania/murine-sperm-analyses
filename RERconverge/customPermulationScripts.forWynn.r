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

