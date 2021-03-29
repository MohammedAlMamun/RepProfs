

###function for sequence alignment and primary processing

PrimaryAnalysis <- function(directory, windowSize){
  
  #function variable defined
  windowSize <- windowSize
  directory <- directory
  
  #set the given directory
  setwd(directory)
  
  #system function calls the command line directives
  system(sprintf(
    
    #directory set up and reading the fastq.gz files
    "cd \"$(pwd)\"
  for sample in `ls *.fastq.gz` 
       do 
       dir=\"$(pwd)\"
       base=$(basename $sample \".fastq.gz\")
       
       #run the alignment
       bowtie2 -x /Users/mohammed/Desktop/Analysis/bowtie2-2.3.5.1-macos-x86_64/indexes/S288C_Ref -U ${dir}/${base}.fastq.gz  -S ${dir}/${base}.sam
       
       #create temporary directory
       temp_dirA=$(mktemp -d)
       
       #obtain sorted bam files
       samtools sort -l 9 -@ 15 -m 1024M  -O bam -o ${dir}/${base}_sorted.bam -T temp_dirA ${dir}/${base}.sam
       
       #index the sorted bam files
       samtools index ${dir}/${base}_sorted.bam
       
       #capture the chromosome size information
       samtools idxstats ${dir}/${base}_sorted.bam | awk 'BEGIN {OFS=\"\\t\"} {if ($2>0) print ($1,$2)}' > ${dir}/${base}_GenomeInfo.txt
       
       #makewindows of appropriate size as defined by the windowSize parameter
       bedtools makewindows -g ${dir}/${base}_GenomeInfo.txt -w %s > ${dir}/${base}_windows.bed
       
       #calculate per nucleotide coverage
       samtools view -h -@ 8 -q 30 -L ${dir}/${base}_windows.bed ${dir}/${base}_sorted.bam | grep -v XS:i: | samtools view -@ 8 -b - | bedtools genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"SRR\",$3}' > ${dir}/${base}.bed
       
       #just creating a directory to save my final bed files
       mkdir dirBeds
       
       #sum the per nucleotide coverage in each bin
       bedtools map -a ${dir}/${base}_windows.bed -b ${dir}/${base}.bed -null 0 -o sum | awk 'BEGIN {OFS=\"\\t\"} {if ($4>=0) print $1,$2,$3,\"SRR\",$4}' > dirBeds/${base}_RpBin.bed
       
       #remove the intermediate files and temporary folders for the sake of memory
       rm ${dir}/${base}.sam
       rm ${dir}/${base}_sorted.bam
       rm ${dir}/${base}_sorted.bam.bai
       rm ${dir}/${base}_GenomeInfo.txt
       rm ${dir}/${base}_windows.bed
       rm ${dir}/${base}.bed
      
       done", windowSize))
  
}

#example run
PrimaryAnalysis("/detail/path/to/folder/containing/your/fastq/files" , windowSize = 100)



###function to make the ratio between chIP and input

MakeRatio <- function(IpCov, InputCov, FACScoefficient){
  
  #require R package Repliscope
  library(Repliscope)
  
  #define function variables
  FACScoefficient <- FACScoefficient
  
  #path link to bed formatted ChIP and Input files
  IpCov <- "path/link/to/the/processed/ChIP/in/bed/format"
  InputCov <- "path/link/to/the/processed/Input/in/bed/format"
  
  #load the bed files
  Ip <- loadBed(IpCov)
  Input <- loadBed(InputCov)
  
  #make ratio between ChIP and Input
  Ratio <- makeRatio(Ip, Input)
  
  #normalise the ratio by the FACS coefficient
  Ratio <- normaliseRatio(Ratio, rFactor = FACScoefficient)
  
  #log10 transformation of the ratio column
  Ratio$ratio <- log10(Ratio$ratio)
  
  #convert any infinite values as 0
  Ratio$ratio[!is.finite(Ratio$ratio)] <- 0
  
  #smooth the ratio using the smoothing function in Repliscope
  Ratio <- smoothRatio(Ratio)
  
  #just removing chrM and remove all unnecessary columns
  F.Ratio <- Ratio[!Ratio$chrom=="chrM", c(1,2,3,9,4)]
  
  #just for ease removing any 0 valued rows
  F.Ratio[F.Ratio==0] <- NA
  F.Ratio <- F.Ratio[complete.cases(F.Ratio), ]
  
  #creating a temporary directory
  to <- tempdir()
  dir.create(to)
  setwd(to)
  
  #saving the ration file as bed 
  write.table(F.Ratio, file="Ratio.bed", quote=FALSE, row.names=FALSE, col.names = FALSE, sep="\t")
  
  #convert the bed to bigwig as deeptool-computematrix considers bigwig input
  system("
    awk '{ print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}' Ratio.bed > Ratio.bedgraph
    sort -k1,1 -k2,2n Ratio.bedgraph > Ratio_sorted.bedgraph
    bedGraphToBigWig Ratio_sorted.bedgraph /Users/mohammed/Desktop/GenomeInfo.txt Ratio.bw
       ")
  
  #renaming the ratio bed as per the given chIP name and saving in a directory called ratio.beds
  file.rename("Ratio.bed", paste0(tools::file_path_sans_ext(basename(IpCov)),".bed"))
  ifelse(!dir.exists(file.path(dirname(IpCov), "ratio.beds")), dir.create(file.path(dirname(IpCov), "ratio.beds")), FALSE)
  file.copy(paste0(tools::file_path_sans_ext(basename(IpCov)),".bed"), paste0(dirname(IpCov), "/ratio.beds") )
  
  #renaming the ratio bigwig as per the given chIP name and saving in a directory called ratio.bigwigs
  file.rename("Ratio.bw", paste0(tools::file_path_sans_ext(basename(IpCov)),".bw"))
  ifelse(!dir.exists(file.path(dirname(IpCov), "ratio.bigwigs")), dir.create(file.path(dirname(IpCov), "ratio.bigwigs")), FALSE)
  file.copy(paste0(tools::file_path_sans_ext(basename(IpCov)),".bw"), paste0(dirname(IpCov), "/ratio.bigwigs") )
  
  #unlink the temporary file
  unlink(to, recursive=TRUE)
  
}
#example run
MakeRatio(IpCov = "/full/path/to/the/processed/chIP/.bed", InputCov = "/full/path/to/the/processed/Input/.bed", FACScoefficient = 1.20)



###plotting the average profiles using computematrix and plotprofile from deeptools
setwd("/directory/info/for/ratio.bigwigs")
dir.create("Matrixes")
dir.create("average.plots")
system("
        computeMatrix reference-point --referencePoint center -b 50000 -a 50000 |
        -R /path/to/bed/formatted/inter-origin/gap/list/category/smallgaps.bed ./../mediumgaps.bed ./../largegaps.bed |
        -S ratio1.bw ratio2.bw ratio3.bw ratio4.bw |
        -o Matrixes/xyz.mat 

        plotProfile --m Matrixes/xyz.mat --numPlotsPerRow 3 --perGroup  --averageType median --plotType se |
        --regionsLabel short mid long --yAxisLabel \"..........\" --yMax 0.40 --yMin 0.00 |
        --samplesLabel xx yy zz ..  --legendLocation upper-right  |
        --o average.plots/xyz.pdf
      ")


###function to plot the local profiles
PlotLocalProfiles <- function(Chromosome, 
                                Ratio_1, Ratio_2, 
                                Ratio_3, Ratio_4, 
                                Start, End, ...){
  #function variables
  #Chromosome <- name of the chromosome to be selected
  #Ratio_files <- processed bed files or dataframe for ratio between chIP and input 
  #Start <- start coordinate
  #End <- end coordinate
  
  #function to select any given chromosomal region from a proceesed bed file or dataframe
  Chr_Choords <- function(ChromDF, Chromosome, Start, End){
    
    #function variables
    #ChromDF <- is the processed dataframe or loaded bedfile
    #Chromosome <- name of the chromosome to be selected
    #Start <- start coordinate
    #End <- end coordinate
    
    #select the given chromosome
    DF <- ChromDF[ChromDF$chrom==Chromosome, ]
    
    #cut from start to end
    DF <- DF[DF$chromStart>=Start & DF$chromEnd<=End, ]
    
    #add a column with the midpoints of standard bed format chromStart and chromEnd
    #to be used for plotting purpose
    DF$Mid <- round((DF$chromStart+DF$chromEnd)/2)
    
    #return the selected region data
    return(DF)
  }
  
  #select origins within that given regions from bed formatted origin list
  Origins <- Chr_Choords(E_Ori, Chromosome, Start, End)
  
  #ratio files are cut and processed for plotting;
  #number ration files are according to experiment & I have four
  Ratio_1 <- Chr_Choords(Ratio_1, Chromosome, Start, End)
  Ratio_2 <- Chr_Choords(Ratio_2, Chromosome, Start, End) 
  Ratio_3 <- Chr_Choords(Ratio_3, Chromosome, Start, End)
  Ratio_4 <- Chr_Choords(Ratio_4, Chromosome, Start, End)
  
  #call plot function from R base and just plot : )
  #I plotted smoothedratio against the midvalue column
  plot(Ratio_1$Mid, Ratio_1$smoothedratio, 
       type="p", col="black", 
       cex=0.5, ylim=c(0,0.4), 
       xlab = " ", ylab = " ",
       main = Title, col.main="black", 
       font.main = 3, cex.main=1.5)
  points(Ratio_2$Mid, Ratio_2$smoothedratio, type="p", col="blue", cex=0.5)
  points(Ratio_3$Mid, Ratio_3$smoothedratio, type="p", col="green", cex=0.5)
  points(Ratio_4$Mid, Ratio_4$smoothedratio, type="p", col="red", cex=0.5)
  
  #put a line at the origin position
  abline(v=c(Origins$Mid))
  
  #add origin names
  text(x = c(Origins$Mid)+((End - Start)*0.05), 
       y = 0, xpd = TRUE, labels = Origins$name,
       col="maroon", cex=0.75, srt = 0, font = 3)
  
}
#example run
#defining plot parameters through par
par(mfrow = c(2,2), oma = c(4,3,2,2) + 0.1, mar = c(3,2,3,2) + 0.1)
#run the function
PlotleLocalProfiles(Chromosome = "chrII", 
                    Ratio_1 = ..., 
                    Ratio_2 = ..., 
                    Ratio_3 = ..., 
                    Ratio_4 = ..., 
                    Start = 400000, 
                    End = 500000)


###############



