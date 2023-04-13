#bash one line sed command for removing loci with no trees from trees file
sed '/ENSMUSP.*AAreplace\t$/d' trees_for_RERconverge.tree > trees_for_RERconverge.noBlanks.tree
