#PURPOSE: Groups "clusters" from Green et al. 2018 scRNAseq dataset together into broader cell type categories
#Separating on this fine a scale might not be meaningful and reduce power (i.e., we don't need ~10 different RS clusters; we can group all the cell type clusters that fall under RS together)

#Spermatogenesis cell types (based on Green et al. 2018 Fig. 2B
gc1<-scan("prot_list.GCcluster1.txt", what=character())
gc2<-scan("prot_list.GCcluster2.txt", what=character())
gc3<-scan("prot_list.GCcluster3.txt", what=character())
gc4<-scan("prot_list.GCcluster4.txt", what=character())
gc5<-scan("prot_list.GCcluster5.txt", what=character())
gc6<-scan("prot_list.GCcluster6.txt", what=character())
gc7<-scan("prot_list.GCcluster7.txt", what=character())
#gc8<-scan("prot_list.GCcluster8.txt", what=character()) #No markers, also kinda looks like junk in their UMAP
gc9<-scan("prot_list.GCcluster9.txt", what=character())
gc10<-scan("prot_list.GCcluster10.txt", what=character())
gc11<-scan("prot_list.GCcluster11.txt", what=character())
gc12<-scan("prot_list.GCcluster12.txt", what=character())

write(gc1, file="prot_list.greenEtal2018.spermatogonia.txt", ncolumns=1)
write(union(gc2, gc3), file="prot_list.greenEtal2018.prelep.txt", ncolumns=1)
sc<-Reduce(union, list(gc4, gc5, gc6, gc7))
write(sc, file="prot_list.greenEtal2018.spermatocytes.txt", ncolumns=1)
st<-Reduce(union, list(gc9, gc10, gc11))
write(st, file="prot_list.greenEtal2018.spermatids.txt", ncolumns=1)
write(gc12, file="prot_list.greenEtal2018.elongating.txt", ncolumns=1)

#Spermatogonia cell types
#spg1<-scan("prot_list.SPGcluster1.txt", what=character())
#spg2<-scan("prot_list.SPGcluster2.txt", what=character())
#spg3<-scan("prot_list.SPGcluster3.txt", what=character())
#spg4<-scan("prot_list.SPGcluster4.txt", what=character())
#spg<-Reduce(union, list(spg1, spg2, spg3, spg4))
#write(spg, file="prot_list.spermatogonia.txt", ncolumns=1)

#Somatic cell types
end<-scan("prot_list.SomaticEndothelial.txt", what=character())
mac<-scan("prot_list.SomaticMacrophage.txt", what=character())
unk<-scan("prot_list.SomaticUnknown.txt", what=character())
lym<-scan("prot_list.SomaticInnateLymphoid.txt", what=character())
myo<-scan("prot_list.SomaticMyoid.txt", what=character())
ley<-scan("prot_list.SomaticLeydig.txt", what=character())
ser<-scan("prot_list.SomaticSertoli.txt", what=character())
somatic<-Reduce(union, list(end, mac, unk, lym, myo, ley, ser))
write(somatic, file="prot_list.greenEtal2018.somatic.txt", ncolumns=1)

print("Done with 04_group_singleCell_clusters.r")
