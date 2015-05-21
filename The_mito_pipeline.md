# Pulling-out-mitogenomes-from-UCE-data
It ain't pretty, but it (probably) works...

1) Create your working directory and copy over the contigs from your trinity-assemlies folder. Obviously change these to your own directories and paths. Once you are inside your trinity-assemblies folder, this is pulling out the Trinity.fasta file within each species folder and renaming it species_name.fasta (based on your folder names)
```
mkdir froggie-mtDNA
dest=/home/a499a400/Kaloula/froggie-mtDNA
#cd to trinity-assemblies folder
for i in `ls`; do cp $i/Trinity.fasta $dest/$i.fasta; done;
```

2) create a blast database of your 'reference' mtDNA. I downloaded blast into my own directory to do this on complab2, but I bet you could do it on the cluster. Basically, I pulled all the Kaloula mtDNA off Genbank that I could find and downloaded it as a fasta file (All_mtDNA_kaloula_fasta), but you could just restrict this to individual species' mitogenomes (or genes e.g. cytb) if you wanted to. Genbank fasta downloads are awful and have line breaks in them, so the grep stage is just to make sure we have no blank lines in them which upsets BLAST.
```
grep -v '^$' All_mtDNA_kaloula.fasta > All_mtDNA_kaloula_no_blank_lines.fasta
/home/a499a400/bin/ncbi-blast-2.2.30+/bin/makeblastdb -in All_mtDNA_kaloula_no_blank_lines.fasta -dbtype nucl
```

3) For each of your fasta files containing the contigs for each sample, we are carrying out a blast search against our mtDNA which is returning the name of the contigs for each sample which map to the mtDNA at 75% % identity (you can make this more or less stringent) 

```
for i in `ls kaloula*`; do /home/a499a400/bin/ncbi-blast-2.2.30+/bin/blastn -db All_mtDNA_kaloula_no_blank_lines.fasta -query $i -perc_identity 75 -outfmt '10 qseqid sseqid qlen slen' | cut -d , -f 1 | uniq > $i.namelist; done;
```

4) You'll need to install seqtk before this step. This is just pulling out the fasta sequence associated with each of the contigs we "hit", so we can use it as a new blast database, to repeat the step again.  We can then get rid of namelist as well, because our info is in oldblast.fasta

```
for i in `ls *.namelist`; do seqname=`echo $i | sed "1s/.fasta.namelist//"` ; ~/bin/seqtk/seqtk subseq $seqname.fasta $i | sed 's/>/>'"$seqname"'+/g' >> $seqname.oldblast.fasta; done;
rm -rf *.namelist
```

5) We are using our new dataset of "sample specific" mtDNA from our actual contigs to do another blast search, just in case we missed anything using the more generic reference last time. Then pulling out the fasta data using seqtk

```
for i in `ls *.oldblast.fasta`; do seqname=`echo $i | sed "1s/.oldblast.fasta//"` ; /home/a499a400/bin/ncbi-blast-2.2.30+/bin/makeblastdb -in $i -dbtype nucl ;  /home/a499a400/bin/ncbi-blast-2.2.30+/bin/blastn -db $i -query $seqname.fasta -perc_identity 95 -outfmt '10 qseqid sseqid qlen slen' | cut -d , -f 1 | uniq > $seqname.namelist; done;
for i in `ls *.namelist`; do seqname=`echo $i | sed "1s/.namelist//"` ; ~/bin/seqtk/seqtk subseq $seqname.fasta $i | sed 's/>/>'"$seqname"'+/g' >> $seqname.newblast.fasta; done;
```

6) Alright - here is where I probably should have got tricky and tried to write some kind of loop (but I was too lazy). If newblast is bigger in size than oldblast for any of the samples (run ls -l to see), delete old blast, rename newblast oldblast, and start again at the rm -rf *.namelist command (last part of #4) until you pull out no more sequences (you can use the code below to do the renaming). If newblast is bigger than oldblast, basically what we are saying is that you pulled out additional contigs on that second blast step (step #5), so we probably should keep going until you pull no more out! You could do this just for the problem children samples, but it doesn't take long just to do it for everybody...

```
rm -rf *.oldblast.*
for i in `ls *.newblast.fasta`; do seqname=`echo $i | sed "1s/.newblast.fasta/.oldblast.fasta/"` ; mv $i $seqname; done ;
```

7) If you don't care about coverage, you are done! This should be all the mtDNA-like contigs for each species, given as species_name.newblast.fasta. If you do care, read on!

8) At this stage, I think it is a good idea to download each of the *.newblast.fasta files and assemble each one against your favorite mitogenome reference. You can do this in the freebie version of Geneious. I just like to check that none of the contigs overlap (in which case, potentially you should think about merging them?).

9) OK, you've done step 8 and merged (or not) your contigs and reuploaded your *.newblast.fasta file (or not - just used the one you've already got there). We are now going to use an R script to buffer each contig with Ns on each side (I've found the assembler we are going to use - bwa - sort of "erodes" the ends of contigs if you don't pad the end of them - probably because we are trying to map PE reads. I'm padding it with 100 Ns, but you can tweak it if you would like!). I also feel like the contig names are super long, so have simplified those a little using this script. Save the following R script as Ning.R and upload it into the same folder as your *.newblast.fasta files. This is going to assume you have previously installed the R package 'stringr'

```
library(stringr)
intable <- readLines("temp")

rows <- length(intable)
to_write <- matrix(NA,ncol=1,nrow=rows)

Ns <- "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

for (j in 1:rows) {
if ((length(grep(">",intable[j])))>0) {
to_write[j] <- unlist(strsplit(intable[j],"\t"))[1]
} else {
to_write[j] <- paste(Ns,intable[j],Ns,sep="")
}
}

write.table(to_write, "tempout",quote=FALSE, col.names=FALSE,row.names=FALSE)
q()
```

10) To actually pad all your contigs in all your files, run the following. This step is slightly untested because I didn't do this in my original foray into the Kaloula mitogenomes (but this is how I would do it now!), so let me know if you run into any problems.

```
for i in *.newblast.fasta; do rm temp*; mv $i temp; Rscript Ning.R; mv tempout $i; done;
```

11) Now we need to fetch our cleaned reads so we can assemble them. First, cd into your cleaned-reads directory, and then do the following (make sure to change your directory!)

```
for i in `ls`; do cp $i/split-adapter-quality-trimmed/$i-READ1.fastq.gz ~/Kaloula/froggie-mtDNA/; cp $i/split-adapter-quality-trimmed/$i-READ2.fastq.gz ~/Kaloula/froggie-mtDNA/; done;
```

12) This takes you through bwa for each of your samples, and through samtools dedupping, and finally, is spitting out both a "straight" pileup (all sites present, no matter how low in coverage) and a "min.10" pileup - only spitting out sites with 10+ coverage. You can tweak this threshold for your own purposes (irl I actually used a threshold of four).

```
for i in `ls *.newblast.fasta` ; do seqname=`echo $i | sed "1s/.newblast.fasta//"` ; ~/bin/bwa/bwa index -a is $i; ~/bin/samtools/bin/samtools faidx $i; ~/bin/bwa/bwa mem $i  $seqname-READ1.fastq.gz $seqname-READ2.fastq.gz > $seqname-alnpe.sam; ~/bin/samtools/bin/samtools view -bS -o $seqname-alnpe.bam $seqname-alnpe.sam; ~/bin/samtools/bin/samtools sort $seqname-alnpe.bam $seqname-alnpe-sorted; ~/bin/samtools/bin/samtools rmdup -s $seqname-alnpe-sorted.bam $seqname-alnpe-sorted-rmdup.bam; ~/bin/samtools/bin/samtools mpileup -f $i $seqname-alnpe-sorted-rmdup.bam > $i.all.pileup; cat $i.all.pileup | awk '{if($4>=10) print $1 " " $2 " " $3 " " $4 " " $5}' > $i.min10.pileup; rm -rf $seqname-aln* ; done;
```

13) If you are interested in calculating average depth etc., the pileup files open pretty nicely in excel and you can average this, otherwise (I think - again, untested) this R script should convert the pileup into fasta (for each sample, every contig will have a separate fasta file. Let me know if this isn't what you want and I can tweak it). Save as converttofasta.R and upload into the directory with all the pileups.

```
intable <- read.table("temp", header=FALSE, stringsAsFactors=FALSE,sep="\t")
maxlength <- dim(intable)[1]

for (j in 1:maxlength) {

commas <- nchar(intable[j,5])-nchar(gsub(",","",intable[j,5]))
dots <- nchar(intable[j,5])-nchar(gsub("\\.","",intable[j,5]))
stars <- nchar(intable[j,5])-nchar(gsub("\\*","",intable[j,5]))
plusses <- nchar(intable[j,5])-nchar(gsub("\\+","",intable[j,5]))
minuses <- nchar(intable[j,5])-nchar(gsub("-","",intable[j,5]))
As <- ((nchar(intable[j,5]))*2)-nchar(gsub("a","",intable[j,5]))-nchar(gsub("A","",intable[j,5]))
Cs <- ((nchar(intable[j,5]))*2)-nchar(gsub("c","",intable[j,5]))-nchar(gsub("C","",intable[j,5]))
Gs <- ((nchar(intable[j,5]))*2)-nchar(gsub("g","",intable[j,5]))-nchar(gsub("G","",intable[j,5]))
Ts <- ((nchar(intable[j,5]))*2)-nchar(gsub("t","",intable[j,5]))-nchar(gsub("T","",intable[j,5]))

if(intable[j,3]=="A") {
As <- As + commas + dots
}
if(intable[j,3]=="C") {
Cs <- Cs + commas + dots
}
if(intable[j,3]=="G") {
Gs <- Gs + commas + dots
}
if(intable[j,3]=="T") {
Ts <- Ts + commas + dots
}

intable[j,6] <- As
intable[j,7] <- Cs
intable[j,8] <- Gs
intable[j,9] <- Ts
intable[j,10] <- stars + plusses + minuses
intable[j,11] <- max(intable[j,6:10])/intable[j,4]

whichmax <- which((intable[j,6:10])==(max(intable[j,6:10])),arr.ind=TRUE)[2]

if(whichmax==1) {
intable[j,12] <- "A"
}
if(whichmax==2) {
intable[j,12] <- "C"
}
if(whichmax==3) {
intable[j,12] <- "G"
}
if(whichmax==4) {
intable[j,12] <- "T"
}
if(whichmax==5) {
intable[j,12] <- "INDEL"
}
}

j <- 1

while (j <= maxlength) {
sequence1 <- intable[j,12]
name <- intable[j,1]
filename <- paste(name,"_out.fas",sep="")
header <- paste(">",name,"_out.fas",sep="")

k <- j + 1

while (k <= maxlength) {
if(intable[j,1]==intable[k,1]) {
if (k <= maxlength) {
sequence1 <- paste(sequence1,intable[k,12],sep="")
k <- k + 1
}
}
if((k <= maxlength) && (!(intable[j,1]==intable[k,1]))) {
break
}
}

out <- rbind(header,sequence1)
write.table(out, filename, quote=FALSE, row.names=F,col.names=F)
j <- j + nchar(sequence1)
}

q()
```

14) To run it (change the regex for the pileup if you want to use the .min10.pileup instead)

```
for i in *all.pileup; do rm temp*; cp $i temp; Rscript converttofasta.R; done;
```

Finally) Just in case our little "buffer with Ns" worked and we extended out contigs, probably worth chucking them through Geneious etc again to make sure none have managed to get into merging territory. Also you'll want to search your contigs for "INDEL". bwa doesn't do very well with indels, so you will probably want to flag these for future investigation (i.e. tweak the reference for these and run them manually).
