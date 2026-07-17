r <- read.delim("all_tips.tsv", stringsAsFactors=FALSE)
# collapse donors that share an identical tip SET (codename vs PD id are the same tree)
sets <- tapply(r$tip, r$donor, function(z) paste(sort(unique(z)), collapse="|"))
key  <- sets[r$donor]
# canonical name per tip-set: prefer a PD-style id
canon <- tapply(names(sets), sets, function(d) { pd <- grep("^PD[0-9]", d, value=TRUE)
                                                 if (length(pd)) pd[1] else d[1] })
r$canon <- canon[key]
r$alias <- ifelse(r$canon==r$donor, "", r$donor)
u <- unique(r[,c("canon","tip")])
n <- table(u$canon)
dep <- tapply(r$median_depth, r$canon, function(z) suppressWarnings(median(z, na.rm=TRUE)))
al  <- tapply(r$alias, r$canon, function(z){ z<-unique(z[z!=""]); paste(z, collapse=",") })
out <- data.frame(donor=names(n), tips=as.integer(n), median_depth=round(dep[names(n)]),
                  alias=al[names(n)], stringsAsFactors=FALSE)
out <- out[order(-out$tips),]
write.table(out, "donors_dedup.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(unique(r[,c("canon","tip")]), "canon_tips.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("distinct donors after dedup: %d  (was %d)\n", nrow(out), length(unique(r$donor))))
