#libraries
library(zoo)
library(seqinr)
library(stringr)
library(bioseq)
library(Biostrings)
library(stringdist)
library(stringi)
library(random)


#import the selected kmers(csv file with manually picked probes, must have a column named 'sequence_with_Ns')
kmers<-read.csv2('C:/Users/nerit/Desktop/new_hits.csv')

#import database (transcriptomes, fasta file with all the transcriptomes)
db<-readDNAStringSet('C:/Users/nerit/Desktop/refmouse.fa', format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

##########################searching algorithm
#function that searches a sequence x for matches in the transcriptomes db we have created, the probes are degenerated, meaning they include N universal binding artificial bases, those match to everything [therefore the argument fixed=FALSE is passed])
#the vmatchpatter function works itteratively checking one query sequence for matches in many target sequences (db)
search_degenerated_probe<-function(x){
  matches<-vmatchPattern(x,db, max.mismatch=5)
}

#function that transforms a regular character string to a DNA string
to_dna_string<-function(x){
  seq<-DNAString(x)
  return(seq)
}
#we get the query(degenerated probe sequences that we will check for matches in the db, we turn those sequences into DNA strings so that the search function [vmatchpattern(pattern,subject,fixed=FALSE)])
query<-sapply(kmers$sequence_with_Ns,FUN=to_dna_string)

#vectorization so that the search function works on all the queries (probes) we have in the manually created file
matches<-sapply(query,FUN=search_degenerated_probe)

#in this for loop we do 3 things: first we remove all the Null element of the "matches" variable (i.e. the non-matches), then create another column in the kmers file, thats the hit column that has all the hits of each query in all the database sequences, separated by ';' and third we introduce another column that counts those hits)
for (i in 1:length(matches)){
  matches[[i]]<-matches[[i]][!sapply(startIndex(matches[[i]]), is.null)] 
  kmers$hits[i]<-paste(matches[[i]]@NAMES,collapse=';')
  kmers$number_of_hits[i]<-length(matches[[i]]@NAMES)
}

#save scv2
write.csv2(kmers,'C:/Users/nerit/Desktop/hits_on_mouse.csv',row.names = FALSE)

