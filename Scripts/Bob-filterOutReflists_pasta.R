## The object with the reflib / Genbank ID indexing is a list called all.taxs.source.otu
## it is larger than the OTU tables oin your phlyoseq objects because it includes all OTUs
## even ones that had no classification from Genbank that were thrown out immediately
## The last column of each data frame in the list is $idsource - here you can see where the OTU was identified from

# View(all.taxs.source.otu)
# View(all.taxs.source.otu$dada.nopool.nochim)
## You can use the OTU code to index and join the data frames.. 

Ref.lib.ids <- lapply(all.taxs.source.otu, function(x) subset(x, idsource == "RefLib"))
# rownames(Ref.lib.ids[[5]])
# View(Ref.lib.ids$dada.pspool.nc.lulu)
## Note that OTU34 needs removal - it was from an unidentified "vine" (i.e. paraphyletic)
## and thus has identification use at all
Ref.lib.ids <- lapply(Ref.lib.ids, function(x) x[!(row.names(x) %in% "OTU34"), ])

# View(Ref.lib.ids)
# View(Ref.lib.ids$dada.pspool.nc.lulu)

## To filter your phyloseq objects 
reflib.otus <- lapply(Ref.lib.ids, function(x) rownames(x))

## To filter your phyloseq objects 
reflib.otus <- lapply(Ref.lib.ids, function(x) rownames(x))

## OK so now the names your Taxa found in the soil that was identified in the reference libraries 
## is in the object reflib.otus - a list of 6 charcater vectors for: dada.no.pool, dada.pooled, dada.psuedo pooled, lulu.no.pool, lulu.pooled, lulu.psuedo.pooled

## So you make your new lists with (for example, on the R1phylo.plants.tax object:
## bob I know you would write a loop because you are a wild and crazy guy...
cutlist1 <- prune_taxa(reflib.otus[[1]], R1phylo.plants.tax[[1]])
cutlist2 <- prune_taxa(reflib.otus[[2]], R1phylo.plants.tax[[2]])
cutlist3 <- prune_taxa(reflib.otus[[3]], R1phylo.plants.tax[[3]])
cutlist4 <- prune_taxa(reflib.otus[[4]], R1phylo.plants.tax[[4]])
cutlist5 <- prune_taxa(reflib.otus[[5]], R1phylo.plants.tax[[5]])
cutlist6 <- prune_taxa(reflib.otus[[6]], R1phylo.plants.tax[[6]])

R1ref.lib.list <- list(cutlist1, cutlist2, cutlist3, cutlist4, cutlist5, cutlist6 )
names(R1ref.lib.list) <- c("R1.dada.nopool.reflib", "R1.dada.pool.reflib" , "R1.dada.pspool.reflib", 
                           "R1.lulu.nopool.reflib", "R1.lulu.pool.reflib", "R1.lulu.pspool.reflib" )

## Getting how many sequences are cut out when only including OTUs
precutsums <- rowSums(otu_table(R1phylo.plants.tax[[6]]))
postcutsums <- rowSums(otu_table(cutlist6))

plot(precutsums, postcutsums, xlab ="All Plants OTUs - library size", ylab ="Only Ref. Lib. OTUs - library size ")
abline(0,1, col="red")

beforeseqs <- lapply(R1phylo.plants.tax, function(x) sum(otu_table(x)))
afterseqs <- lapply(R1ref.lib.list, function(x) sum(otu_table(x)))
unlist(afterseqs)/unlist(beforeseqs)



## I am not sure that it is relevant to have an "other" column with everything else collapsed other than to say that x% of reads from plants matched to the reference
## library and for this all you need to do is compare the vectors from rowSums for the pre-reflib filtered otu_table() and the post-reflib filtered otu_table()

## Now to link the OTU names to the POTU table that cesc has, just access the object:

## to get the link between otu names and POTU (cesc's code for the reference libraries, load this list)
otu.potu.link <- readRDS("/Users/glennd/Documents/GitHub/legenDNAry/Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_v1.rds")

colnames(otu_table(R1ref.lib.list[[5]]))
otu.potu.link[[5]]$OTU

View(R1ref.lib.list)


###
R1ref.lib.list[[6]]
