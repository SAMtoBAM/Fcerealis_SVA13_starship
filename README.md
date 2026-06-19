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
#SVA13 = GCA054574715 
#F2344 = GCA054553065

for longread in GCA054574715 GCA054553065
do
ls Fcerealis_public_assemblies/${longread}.fa
done > assemblies_list.LR.txt 

##run stargraph but now using a reduced assemblies list
stargraph.sh -a assemblies_list.LR.txt -t 16 -p Fcerealis -o stargraph_output -r starfish_wrapper_output/starfish.filt.SRGs_combined.gff -e  starfish_wrapper_output/starfish.elements.ann.FILTERED.feat

```


### Step 1d: Run cargobay to find Horizontal gene transfer

```bash



```
