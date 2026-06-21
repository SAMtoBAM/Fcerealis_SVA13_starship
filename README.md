# Fcerealis_SVA13_starship
Repository for the analysis of Starships in the Fusarium cerealis strain SVA13 genome

Using the assembly of SVA13 made publicly available for ease of reproducibility <br/>
SVA13 = GCA_054574715.1 <br/>

Other long-read assembly is F23-4.4 (GCA_054553065.1)

## Step 0: set up environment

```bash
##going to be mainly using the stargraph tool conda environment

conda create -n starg samtobam::stargraph
conda activate starg
##will just add some R packages for the phylogeny too
conda install bioconda::bioconductor-ggtree conda-forge::r-ape conda-forge::r-phangorn conda-forge::r-ggplot2

##create and move into directory for analysis
mkdir Fcerealis_starships && cd Fcerealis_starships

```


## Step 1: Identify _Starships_ and SLRs in all public Fusarium cerealis assemblies; including 2 long-read assemblies for SVA13 and F-23-4.4 (F2344)
### Step 1a: Download all assemblies and rename

```bash

mkdir Fcerealis_public_assemblies

##download the assembly and protein dataset then rename it etc
datasets download genome taxon "Fusarium cerealis"
unzip ncbi_dataset.zip
rm ncbi_dataset.zip
##rename them as just the ncbi GCA assession and rename the contig headers with the genome name 
ls ncbi_dataset/data/ | grep -v json | while read genome
do
    genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
    cat ncbi_dataset/data/$genome/$genome*.fna | sed "s/>/>${genome2}_/g" | awk -F " " '{print $1}' | awk '{if($0 ~ ">") {print} else {print toupper($0)}}' > Fcerealis_public_assemblies/$genome2.fa
done
##clean up
rm -r ncbi_dataset/ md5sum.txt README.md 
```

### Step 1b: run starfish_wrapper (mainly just to get the de-novo Starship gene annotations


```bash
##get list of assemblies
ls Fcerealis_public_assemblies/* > assemblies_list.txt

##launch stargraph's starfish_wrapper
starfish_wrapper.sh -a assemblies_list.txt -t 16 -o starfish_wrapper_output

##no starships detected so create a empty file for stargraph mimicing the output (just called it what I might call it if there were starships)
touch starfish_wrapper_output/starfish.elements.ann.FILTERED.feat

```

### Step 1c: Run stargraph on long-read assemblies SVA13 and F-23-4.4

```bash
##now only want to run stargraph with the long-read assemblies
SVA13="GCA054574715" 
F2344="GCA054553065"

for longread in ${SVA13} ${F2344
do
ls Fcerealis_public_assemblies/${longread}.fa
done > assemblies_list.LR.txt 

##run stargraph but now using a reduced assemblies list
stargraph.sh -a assemblies_list.LR.txt -t 16 -p Fcerealis -o stargraph_output -r starfish_wrapper_output/starfish.filt.SRGs_combined.gff -e  starfish_wrapper_output/starfish.elements.ann.FILTERED.feat

```
Raw output for a Starship in SVA13 <br/>
Clearly the edges is a little off but we can refine that later. Captin edge of predicted element extends over aligned region, and other end extends over ~200bp of telomeric sequences. Also as we will see below the 20kb prior to the telomeres is likely not Starship Material neither.

<p>
<img src="https://github.com/SAMtoBAM/Fcerealis_SVA13_starship/blob/main/images/cluster2.svg" width=100%>
</p>


## Step 2: HGT evidence 
### Step 2a: Run cargobay to find Horizontal gene transfer (currently not working as database is down)

```bash
##first need to create a metadata file that will tell cargobay which species each genome is and therefore not consider it HGT
for longread in GCA054574715 GCA054553065
do
echo "${longread};Fusarium cerealis" | tr ';' '\t'
done > cargobay_metadata.tsv

cargobay.sh -t 2 -p Fcerealis -o cargobay_output -a stargraph_output/Fcerealis.assemblies.fa.gz -e stargraph_output/Fcerealis.starships_SLRs.fa -b stargraph_output/Fcerealis.starships_SLRs.bed -m cargobay_metadata.tsv -g starfish_wrapper_output/starfish.filt.SRGs_combined.gff

```


### Step 2b: Manually identify genomes with the SVA13 starship and generate alignments with Starship gene predictions

```bash
##just want to extract the starship region pus some flank (can only go further downstream)
samtools faidx stargraph_output/Fcerealis.assemblies.fa.gz $( grep $SVA13 stargraph_output/Fcerealis.starships_SLRs.bed | awk '{print $1":"$2"-"$3+20000}' ) > SVA13.GCA054574715_SLR2.20kb_flank.fa
```

Based on BLASTing NCBI we can find this element in 3 genomes of other species: <br/>
Fusarium culmorum = GCA_052570865.1 <br/>
Fusarium sp. = GCA_022627095.1 <br/>
Fusarium equiseti = GCA_019055085.1

```bash

candidates="GCA_052570865.1 GCA_022627095.1 GCA_019055085.1"

mkdir HGT_candidates_genomes

##download the assembly and protein dataset then rename it etc
datasets download genome accession ${candidates}
unzip ncbi_dataset.zip
rm ncbi_dataset.zip
##rename them as just the ncbi GCA assession and rename the contig headers with the genome name 
ls ncbi_dataset/data/ | grep -v json | while read genome
do
    genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
    cat ncbi_dataset/data/$genome/$genome*.fna | sed "s/>/>${genome2}_/g" | awk -F " " '{print $1}' | awk '{if($0 ~ ">") {print} else {print toupper($0)}}' > HGT_candidates_genomes/$genome2.fa
done
##clean up
rm -r ncbi_dataset/ md5sum.txt README.md 

##combine them
cat HGT_candidates_genomes/*.fa > HGT_candidates_genomes.fa

##BLAST them with the starship plus flank and extract the aligning contigs (with a minimum of 5kb alignment)
##this gives us 1 contig per genome (as expected); and then we can aggregate the alignments to identify a large contiguous aligned region
##and place it in a tsv
echo "contig;start;end;size" | tr ';' '\t' > HGT_candidates_genomes.contigs.bed
blastn -query SVA13.GCA054574715_SLR2.20kb_flank.fa -subject HGT_candidates_genomes.fa -outfmt 6  |
awk '
$4 >= 5000 { keep[$2]=1 }
{
    lines[NR]=$0
}
END{
    for(i=1;i<=NR;i++){
        split(lines[i],f,"\t")
        if(keep[f[2]])
            print lines[i]
    }
}' |
awk '
{
    s=($9<$10)?$9:$10
    e=($9>$10)?$9:$10
    print $2"\t"s"\t"e
}' |
sort -k1,1 -k2,2n |
awk -v gap=500 -v minlen=10000 '
NR==1{
    chr=$1; start=$2; end=$3
    next
}
{
    if($1==chr && $2-end < gap){
        if($3>end) end=$3
    } else {
        if(end-start+1 >= minlen)
            print chr"\t"start"\t"end"\t"end-start+1
        chr=$1; start=$2; end=$3
    }
}
END{
    if(end-start+1 >= minlen)
        print chr"\t"start"\t"end"\t"end-start+1
}' >> HGT_candidates_genomes.contigs.bed

##extract those entire contigs to a file for alignment
##add the contig in SVA13 first for order
samtools faidx stargraph_output/Fcerealis.assemblies.fa.gz $( grep $SVA13 stargraph_output/Fcerealis.starships_SLRs.bed | awk '{print $1}' ) > HGT_candidates_genomes.contigs.fa
##then the others
seqkit grep -f <(cut -f1 HGT_candidates_genomes.contigs.bed) HGT_candidates_genomes.fa >> HGT_candidates_genomes.contigs.fa
##this file will be used for all vs all alignments

##run starfish on the hGT candidate genomes to get the Starship related genes
ls HGT_candidates_genomes/*.fa > HGT_assemblies_list.txt
starfish_wrapper.sh -a HGT_assemblies_list.txt -t 16 -o HGT_candidates_starfish_wrapper_output
##will use this and the SVA13 starfish run to extract the coordinates of the captain and any other starship related gene that may be present on these contigs

##FIRST PLOTTING FILE (bed file with the region where alignments will be based on BLAST and some flanks, will show the lines representing the contig)
##now add 30kb flank to the edges..if possible
flank=30000

##need to index the contigs so I get easily get the contig sizes, used to calculate max length possible for flank
samtools faidx HGT_candidates_genomes.fa

echo "contig;start;end" | tr ';' '\t' > HGT_candidates_genomes.contigs.flank.bed
##add the SVA13 region
grep $SVA13 stargraph_output/Fcerealis.starships_SLRs.bed | cut -f1-3 | while read -r contig start end
do
flankstart=$( echo $start | awk -v flank="$flank" '{if(($1-flank) < 0) {print "0"} else {print $1-flank}}' )
flankend=$( grep $contig stargraph_output/Fcerealis.assemblies.fa.gz.fai | cut -f2 | awk -v end="$end" -v flank="$flank" '{if((end+flank) > $1){print $1} else {print end+flank}}' )
echo "${contig};${flankstart};${flankend}" | tr ';' '\t'
done >> HGT_candidates_genomes.contigs.flank.bed
##add the candidates
tail -n+2 HGT_candidates_genomes.contigs.bed | cut -f1-3 | while read -r contig start end
do
flankstart=$( echo $start | awk -v flank="$flank" '{if(($1-flank) < 0) {print "0"} else {print $1-flank}}' )
flankend=$( grep $contig HGT_candidates_genomes.fa.fai | cut -f2 | awk -v end="$end" -v flank="$flank" '{if((end+flank) > $1){print $1} else {print end+flank}}' )
echo "${contig};${flankstart};${flankend}" | tr ';' '\t'
done >> HGT_candidates_genomes.contigs.flank.bed

##SECOND PLOTTING FILE (alignment all vs all using nucmer converted to a paf file)
##usually use a minmatch of 100 but removing it due to lower identity and telomeric sequences
nucmer -t 16 --maxmatch --delta HGT_candidates_genomes.contigs.nucmer.delta HGT_candidates_genomes.contigs.fa HGT_candidates_genomes.contigs.fa
paftools.js delta2paf HGT_candidates_genomes.contigs.nucmer.delta > HGT_candidates_genomes.contigs.nucmer.paf


##THIRD PLOTTING FILE (coordinates of starship related genes in the aligned contigs)
##first header
echo "contig;start;end;sense;gene;label" | tr ';' '\t' > HGT_candidates_genomes.contigs.genes.bed
##predicted starship related genes for SVA13
cat starfish_wrapper_output/geneFinder_*/starfish.filt.gff | grep "$SVA13"  | awk '{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' | tr ';' '\t' | sed 's/Name=//g' | cut -f1,2,3,4,6 | awk '{if($5 ~ "_duf3723") {print $0"\tDUF3723"} else if($5 ~ "_tyr") { print $0"\ttyrR"} else if($5 ~ "_myb") { print $0"\tMYB"}}' >> HGT_candidates_genomes.contigs.genes.bed
##predicted starship related genes for HGT candidates
cat HGT_candidates_starfish_wrapper_output/geneFinder_*/starfish.filt.gff | awk '{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' | tr ';' '\t' | sed 's/Name=//g' | cut -f1,2,3,4,6 | awk '{if($5 ~ "_duf3723") {print $0"\tDUF3723"} else if($5 ~ "_tyr") { print $0"\ttyrR"} else if($5 ~ "_myb") { print $0"\tMYB"}}' >> HGT_candidates_genomes.contigs.genes.bed

##manually adding the identified Haemolysin III (AdipoR/Haemolysin-III-related (IPR004254)) and RHOD (Rhodanese-like domain (IPR001763)) proteins (based on ORF predictions and interpro searching)
##likely involed in antagnosism and cyanide detoxification respectively and both associated with effectors
echo "GCA054574715_JBJHEB010000022.1;35011;35700;-;HlyIII;HlyIII" | tr ';' '\t' >> HGT_candidates_genomes.contigs.genes.bed
echo "GCA054574715_JBJHEB010000022.1;25998;26732;+;RHOD;RHOD" | tr ';' '\t' >> HGT_candidates_genomes.contigs.genes.bed
##manually identified proteins with singalP domains
echo "GCA054574715_JBJHEB010000022.1;37795;38166;-;putative-effector;putative-effector" | tr ';' '\t' >> HGT_candidates_genomes.contigs.genes.bed
echo "GCA054574715_JBJHEB010000022.1;25654;26037;-;putative-effector;putative-effector" | tr ';' '\t' >> HGT_candidates_genomes.contigs.genes.bed
echo "GCA054574715_JBJHEB010000022.1;21708;22019;+;putative-effector;putative-effector" | tr ';' '\t' >> HGT_candidates_genomes.contigs.genes.bed
##and a Patatin-like phospholipase protein (as found often in Starships)
echo "GCA054574715_JBJHEB010000022.1;21040;21594;+;PLP;PLP" | tr ';' '\t' >> HGT_candidates_genomes.contigs.genes.bed	
	
##NOT CURRENTLY USING
##gene sequences from anntoation
#grep CONTIG GENOME.gff | grep CDS | awk -F ";" '{print $1}' | sed 's/ID=cds-//g' | awk '{print $1"\t"$4"\t"$5"\t"$7"\t"$9"\tNA"}' >> HGT_candidates_genomes.contigs.genes.bed

##FORTH AND FINAL PLOTTING FILE (bed file showing the position of the predicted element
##the edges were manually modified using the alignments with HGT candidates
#SVA13 starship = GCA054574715_JBJHEB010000022.1:19376-44853
##currently naming that starship what is what identified as from stargraph
echo "contig;start;end;starship" | tr ';' '\t' > HGT_candidates_genomes.contigs.starship.bed
echo "GCA054574715_JBJHEB010000022.1;19376;44853;GCA054574715_SLR2" | tr ';' '\t' >> HGT_candidates_genomes.contigs.starship.bed


```
### Step 2c: Plotting alignments

```R

##can use R inside the stargraph conda environment to have all the libraries required, however it is just gggenomes, IRanges and ggnewscale
library(IRanges)
library(gggenomes)
library(ggnewscale)

##first feature
##region to be plotted
bed=read.csv("HGT_candidates_genomes.contigs.flank.bed", sep='\t', header=T)
##just modify some headers for downstream handling
##make another header, copy of contig, but called seq_id
bed$seq_id = bed$contig
##add another column, bin_id, which will be used by gggenomes to split up each element onto its own line (defaults to seq_id for clustering per row if not present)
bed$bin_id = 1:nrow(bed) 
##also length using the end position
#bed$length = bed$end
bed$length = bed$end - bed$start
##add a column with label to be used in the plot (showing the actual region being aligned)
bed$label= paste(bed$contig,":",bed$start,"-",bed$end, sep = "")
##add another column that is just the genome accession
bed$genome <- sub("_.*", "", bed$contig)

##second feature
##a bed file of just the Starship-like region coordinates i.e. the above bed file without the flanking regions
SLRbed=read.csv("HGT_candidates_genomes.contigs.starship.bed", sep='\t', header=T)
SLRbed$seq_id = SLRbed$contig
SLRbed$length = SLRbed$end-SLRbed$start

##third feature
##a bed file with the genes annotated within the Starship-like regions (coordinates modified due to the flanking regions added)
genes=read.csv("HGT_candidates_genomes.contigs.genes.bed", sep='\t', header=T)
genes$seq_id = genes$contig
genes$length= genes$end-genes$start
genes$strand=genes$sense

##fourth feature
##the nucmer all-v-all alignment converted to paf format
links=read_links("HGT_candidates_genomes.contigs.nucmer.paf")


##the actual plot
##reordered for clarity and no filtering on alignments
gggenomes(genes=genes, seqs=bed, feat=SLRbed, links=subset(links, seq_id != seq_id2), adjacent_only = T) %>%
    gggenomes::flip(1) %>%
    gggenomes::sync() %>%
    gggenomes::pick(4,1,3,2) +
    geom_link(aes(fill=((map_match/map_length)*100)) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
    scale_fill_gradientn(colours=c("grey100","grey75", "grey50"), name ="Identity (%)", labels=c(70,80,90,100), breaks=c(70,80,90,100), limits = c(70, 100))+
    new_scale_fill()+
    geom_seq(linewidth = 0.5)+
    geom_feat(color="red", alpha=.6, linewidth=3)+
    geom_gene(aes(fill=label), stroke=0.1, colour="black", shape = 3)+
    geom_seq_label(aes(label=label))+
    geom_seq_label(aes(label=genome), nudge_y = -.25)+
    geom_feat_label(aes(label=starship), nudge_y = -.2, siz=4, angle = 0, fontface = "italic")+
    geom_gene_tag(aes(label=label), size = 2, nudge_y=0.1, check_overlap = FALSE)+
    scale_fill_manual(values = c("red","blue","lightblue","orange", "orange", "orange", "grey"), breaks=c("tyrR","MYB", "DUF3723", "HlyIII", "RHOD","putative-effector", "PLP"), name = NULL)+
    theme(legend.position="top", legend.box = "horizontal")

```
Best candidate as below is the F. culmorum assembly GCA_052570865.1, due to large unaligned flanks and higher identity
<p>
<img src="https://github.com/SAMtoBAM/Fcerealis_SVA13_starship/blob/main/images/SVA13_SLR2_HGT_alignment.svg" width=100%>
</p>



## Step 2e: Phylogeny with HGT candidates

Can generate a quick k-mer based phylogeny using all the F. cerealis genomes, the HGT candidate genomes and some other reference Fusarium species assemblies <br/>
These other reference genomes will include several well known species and several from the sambucinum complex <br/>
Reference genomes used: <br/>
F. oxysporum Fo47 (GCA_013085055.1)  <br/>
F. oxysporum Fo5176 (GCA_030345115.2) <br/>
F. culmorum Class2-1B (GCA_016952355.1)  <br/>
F. culmorum PV (GCA_003033665.1) <br/>
F. poae DAOMC252244 (GCA_019609905.1)  <br/>
F. sambucinum potato_lamoka (GCA_050947815.1)  <br/>
F. graminearum PH-1/NRRL31084 (GCA_000240135.3)  <br/>
F. pseudograminearum CS3096 (GCA_000303195.2) <br/>
F. verticillioides 7600 (GCA_000149555.1) <br/>
F. asiaticum KCTC16664 (GCA_025258505.1) <br/>
F. vorosii RN1 (GCA_037179535.1) <br/>
F. boothii CBS316.73 (GCA_017656985.1) <br/>
F. equiseti S.F-5 (GCA_052857265.1) <br/>
F. venenatum MPI-CAGE-CH-0201 (GCA_020744135.1) <br/>
F. sporotrichioides S17/16 (GCA_019054645.1) <br/>
F. tricinctum MsR-QD66 (GCA_050859235.1) <br/>
F. chlamydosporum IraGTOF6 (GCA_047716405.1) <br/>
F. austroamericanum CBS110246 (GCA_017657035.1) <br/>

```bash

candidates="GCA_013085055.1 GCA_016952355.1 GCA_003033665.1 GCA_019609905.1 GCA_050947815.1 GCA_000240135.3 GCA_000303195.2 GCA_000149555.1 GCA_025258505.1 GCA_037179535.1 GCA_017656985.1 GCA_052857265.1 GCA_020744135.1 GCA_019054645.1 GCA_050859235.1 GCA_047716405.1 GCA_017657035.1"

mkdir assemblies_for_phylogeny

##download the assembly and protein dataset then rename it etc
datasets download genome accession ${candidates}
unzip ncbi_dataset.zip
rm ncbi_dataset.zip
##rename them as just the ncbi GCA assession and rename the contig headers with the genome name 
ls ncbi_dataset/data/ | grep -v json | while read genome
do
    genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
    cat ncbi_dataset/data/$genome/$genome*.fna | sed "s/>/>${genome2}_/g" | awk -F " " '{print $1}' | awk '{if($0 ~ ">") {print} else {print toupper($0)}}' > assemblies_for_phylogeny/$genome2.fa
done
##clean up
rm -r ncbi_dataset/ md5sum.txt README.md 


##gernerate k-mer signatures for each assembly
mkdir signatures

for genome in $( ls */*.fa | grep -v stargraph_output ); do
    sourmash sketch dna \
        "$genome" \
        -p k=21,scaled=1000 \
        --name $(basename ${genome%.*}) \
        -o signatures/$(basename ${genome%.*}).sig
done

##compare the k-mer signatures to generate a similarity matrix
sourmash compare signatures/*.sig -k 21 --distance-matrix --csv kmer_distance.mat

```

Now plot the distance matrix as a midrooted NJ tree in R with a red tip for SVA13 and blue tips for the HGT candidates

```R

library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)

m <- read.csv(
    "kmer_distance.mat",
    header = TRUE,
    check.names = FALSE
)

##get the header and add it to the first column 
taxa <- colnames(m)
m <- as.matrix(m)
rownames(m) <- taxa
colnames(m) <- taxa

# ensure numeric (otherwise get some errors with 0 values)
mode(m) <- "numeric"

# distances already
d <- as.dist(m)

# build tree
tree <- nj(d)

##midroot the tree
tree <- midpoint(tree)

##set a list of the genomes with putative HGT so we can highlight the tips in the tree
SVA13=c("GCA054574715")
HGT=c("GCA052570865","GCA022627095","GCA019055085")

##plot with ggtree
ggtree(tree, size=0.2) +
    guides(fill = "none")+
    geom_tiplab(size=12, as_ylab = TRUE, align = T)+
    scale_fill_viridis_b()+
    geom_rootedge(rootedge = 0.01, linewidth=0.2)+
    geom_treescale(x=0.05, y=12, width=0.05, color='black')+
    geom_tippoint(aes(subset = label %in% SVA13), color = "red", size = 3)+
    geom_tippoint(aes(subset = label %in% HGT), color = "blue", size = 3 )

```
We see that the best aligned clusters with other F. culmorum assemblies <br/>
And the two others are very far away, closest to F. equiseti

<p>
<img src="https://github.com/SAMtoBAM/Fcerealis_SVA13_starship/blob/main/images/SVA13_SLR2_HGT_phylogeny.svg" width=100%>
</p>

Notably this element cannot be comfirmed in F. culmorum using the other culmorums as it appears to be inserted in a region not present in other assembies or the same <br/>
Might be able to add GCA_016952355.1 chromosome2/CP064748 around 7025000 kb

Also likely DR = TTACAG <br/>
and likely TIR = ATAAACCTT <br/>


