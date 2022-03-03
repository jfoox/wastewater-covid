### set up
rootdir=$1
projectname=$2
projectid=$3
mkdir ${rootdir}/${projectname}

# ---------------------------------------------------------------
# this block of code is specific to our internal basespace setup;
# assumes downloading from basespace and merging 4 lanes of
# nextseq output. ultimately all you need is a "lanes_merged"
# folder where every sample is included as:
# ${samplename}_merged.R1.fastq.gz
# ${samplename}_merged.R2.fastq.gz
cd ${rootdir}
basespace download project -i $projectid -o ${projectname}

mkdir ${rootdir}/${projectname}/data_merged
cd    ${rootdir}/${projectname}

# merge lanes
for i in */*fastq.gz; do if [ ! -f $(basename $i) ]; then mv $i .; fi; done 
rm -r *_ds*
parallel -j 20 "cat {}*R1*fastq.gz > data_merged/{}.merged_R1.fastq.gz" ::: $(ls *gz | sed -r "s/_S.+//g" | sort -u)
parallel -j 20 "cat {}*R2*fastq.gz > data_merged/{}.merged_R2.fastq.gz" ::: $(ls *gz | sed -r "s/_S.+//g" | sort -u)
# ---------------------------------------------------------------

numsamples=$(ls ${rootdir}/${projectname}/data_merged/*R1.fastq.gz | wc -l)

# create analysis folder for both
mkdir ${rootdir}/${projectname}
mkdir ${rootdir}/${projectname}/kraken
cd    ${rootdir}/${projectname}/kraken
for i in $(ls ${rootdir}/${projectname}/lanes_merged/*R1*); do echo $i; done | sort -V > ${projectname}_R1s.txt

### run kraken
cd ${rootdir}/${projectname}/kraken
mkdir logs
mkdir outs
# on submission node
sbatch --array 1-$numsamples COVID_kraken.slurm ${projectname}_R1s.txt ${rootdir}/${projectname}/kraken/outs



### gather and visualize kraken outputs
mkdir ${rootdir}/${projectname}/tables
cd    ${rootdir}/${projectname}/tables
echo 'sample,total,unclassified,bacteria,archaea,fungi,human,viruses,sars_cov2' > kraken__taxonomy.csv
for i in ${rootdir}/${projectname}/kraken/outs/*.report.gz; do
  sample=$(basename $i | cut -d'.' -f1)
  reads_classified=$(zgrep 'root$' $i | cut -f2)
  reads_unclassified=$(zgrep -P '\tunclassified' $i | cut -f2)
  reads_total=$(( $reads_classified + $reads_unclassified ))
  reads_bacteria=$(zgrep -P '\tD\t' $i | grep Bacteria | cut -f2)
  reads_archaea=$(zgrep -P '\tD\t' $i | grep Archaea | cut -f2)
  reads_fungi=$(zgrep -P '\tK\t' $i | grep Fungi | cut -f2)
  reads_human=$(zgrep 'Homo sapiens' $i | cut -f2)
  reads_viruses=$(zgrep -P '\tD\t' $i | grep Viruses | cut -f2)
  reads_SARSCoV2=$(zgrep 'Severe acute respiratory syndrome coronavirus 2' $i | cut -f2)
  echo $sample,$reads_total,$reads_unclassified,$reads_bacteria,$reads_archaea,$reads_fungi,$reads_human,$reads_viruses,$reads_SARSCoV2 \
    >> kraken__taxonomy.csv
done

### kraken --> fastqs (sars only)
mkdir ${rootdir}/${projectname}/kraken2fastq
cd    ${rootdir}/${projectname}/kraken2fastq
mkdir logs
mkdir outs
for i in ${rootdir}/${projectname}/kraken/outs/*output.gz; do fq=$(ls ${rootdir}/${projectname}/lanes_merged/*R1* | grep $(basename $i | cut -d'.' -f1)); echo ${fq}@${i}; done | sort -V > fq_and_krakenoutputs.txt
# on submission node
sbatch --array 1-$numsamples COVID_kraken-to-fastq-SARSonly.slurm fq_and_krakenoutputs.txt ${rootdir}/${projectname}/kraken2fastq/outs

### alignment
mkdir ${rootdir}/${projectname}/bwa
cd    ${rootdir}/${projectname}/bwa
function sars_bwa {
  ref_path=/athena/masonlab/scratch/projects/metagenomics/covid/reference/GCF_009858895.2_ASM985889v3_genomic.fna
  R1=$1
  R2=${R1/.R1/.R2}
  sample=$(basename $R1 | cut -d'.' -f1)
  readgroup=@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA
  bwa mem -t 4 -R $readgroup $ref_path $R1 $R2 | samtools view -b - | samtools sort -o ${sample}.bam -
}; export -f sars_bwa
ls ${rootdir}/${projectname}/kraken2fastq/outs/*R1*.fastq.gz \
  | env_parallel -j 4 -v --lb "sars_bwa {}"


### primer trim (assuming ARTIC v3)
cd ${rootdir}/${projectname}/bwa
function ivar_primertrim {
  bam=$1
  sample=$(echo $bam | cut -d'.' -f1)
  bed=artic-nCov-2019.v3.primer.bed
  ivar trim -i $bam -b $bed -e -q 15 -p ${sample}.ivar > ${sample}.ivartrim.log 2> /dev/null
  sambamba sort ${sample}.ivar.bam 2> /dev/null
  mv ${sample}.ivar.sorted.bam ${sample}.ivar.bam 
  mv ${sample}.ivar.sorted.bam.bai ${sample}.ivar.bam.bai
}; export -f ivar_primertrim
ls *.bam | env_parallel --lb -v -j 8 "ivar_primertrim {}"


### mosdepth for amplicon coverage
cd ${rootdir}/${projectname}/bwa
parallel -j 4 "sambamba index {}" ::: *.bam
# v3 samples
parallel -j 4 "mosdepth --threads 3 --by artic-nCov-2019.v3.insert.proper.renamedchr.bed -n {.} {}" ::: *.ivar.bam
mkdir ../mosdepth
mv *mosdepth* ../mosdepth/
mv *regions* ../mosdepth/
cd ../mosdepth
# create table
for i in *gz; do sample=$(echo $i | cut -d'.' -f1); zcat $i | sed -e "s/NC_045512.2/$sample/g"; done | sed -r "1s/^/sample\tstart\tend\tamplicon\tcoverage\n/" > ../tables/mosdepth_allsamples.tsv



### get genome coverages to see if non-mutations are covered or not
cd ${rootdir}/${projectname}/bwa
parallel --lb -v -j 8 "bedtools genomecov -d -ibam {} > {.}.genomecov" ::: *.ivar.bam

# get >=75% samples
cd ${rootdir}/${projectname}/mosdepth
for i in *.bed.gz; do sample=$(echo $i | cut -d'.' -f1); num10x=$(zcat $i | awk '$5 >= 10' | wc -l); if [ $num10x -gt 73 ]; then echo $sample; fi; done > ../samples_genomecov_min75pct.txt

### Call variants: IVAR + LoFreq and take union of calls
mkdir ${rootdir}/${projectname}/varcall
cd    ${rootdir}/${projectname}/varcall

# ivar
function run_ivar {
  bam=$1
  ref=GCF_009858895.2_ASM985889v3_genomic.fna
  sample=$(basename $bam | cut -d'.' -f1)
  
  if $(grep -q $sample ../samples_genomecov_min75pct.txt); then
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference $ref $bam | ivar variants -p ivar_$sample -q 10 -t 0.1 -m 20 -r $ref
    ivar_variants_to_vcf.py ivar_${sample}.tsv ivar_${sample}.vcf
  else echo "$sample not covered well enough"
  fi
}; export -f run_ivar
ls ../bwa/*.ivar.bam | env_parallel --lb -v -j 8 "run_ivar {}"

# lofreq
function run_lofreq {
  bam=$1
  ref=GCF_009858895.2_ASM985889v3_genomic.fna
  sample=$(basename $bam | cut -d'.' -f1)
  
  if $(grep -q $sample ../samples_genomecov_min75pct.txt); then 
    lofreq call --call-indels -f $ref -o lofreq_${sample}.vcf --verbose $bam
  else echo "$sample not covered well enough"
  fi  
}; export -f run_lofreq
ls ../bwa/*.ivar.bam | env_parallel --lb -v -j 20 "run_lofreq {}"

# combine ivar and lofreq VCF info (for importing VAF/DP into R)
for i in ivar*vcf; do
  sample=$(echo $i | cut -d'_' -f2- | cut -d'.' -f1)
  python ../../combine_signal_from_ivar_and_lofreq.py $sample
done


### VEP to annotate variants, on both ivar and lofreq outputs
conda activate vep
for vcf in *.vcf; do
    sample=$(echo $vcf | cut -d'.' -f1)
    vep --verbose --offline --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --mane --protein --species sars_cov_2 --symbol --transcript_version --tsl --input_file $vcf --output_file ${sample}_vep_sarscov2.txt
    python ../../parse_vep_output.py ${sample}
done; rm *warnings.txt; rm *summary.html; rm *_vep_sarscov2.txt
for i in *_vep_sarscov2_reportinput.txt; do mv -i $i ${i/_vep_sarscov2_reportinput.txt/.vepreport.txt}; done

# combine and uniquify ivar and lofreq outputs
for i in ivar*vcf; do
  sample=$(echo $i | cut -d'_' -f2- | cut -d'.' -f1)
  cat *${sample}.vepreport.txt | sort -uV > ivarAndLoFreq_VEP_${sample}.csv
done

# fix headers
for i in ivarAndLoFreq_VEP_*.csv; do cat <(echo "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra,SYMBOL") <(tail -n+2 $i) | grep -v '##' > ${i}_2; done
rm ivarAndLoFreq_VEP_*.csv
for i in *csv_2; do mv -i $i ${i/csv_2/csv}; done


### Freyja for variant demixing
mkdir ${rootdir}/${projectname}/freyja
cd    ${rootdir}/${projectname}/freyja
ref=GCF_009858895.2_ASM985889v3_genomic.fna
for i in ../bwa/*.ivar.bam; do
  sample=$(basename $i | cut -d'.' -f1)
  
  if $(grep -q $sample ../samples_genomecov_min75pct.txt); then
    tsv=../varcall/ivar_${sample}.tsv
    echo [`date`] $i
    freyja variants $i --variants $tsv --depths ${sample}.depths --ref $ref
    if [[ $(wc -l $tsv) > 1 ]]; then 
      freyja demix $tsv ${sample}.depths --output ${sample}.out
    fi
  else echo "$sample not covered well enough"
  fi
done
# make table
for i in *out; do sample=$(echo $i | cut -d'.' -f1); cat $i | grep summarized | cut -f2 | sed -r "s/, \(/\n/g" | sed -e s"/^/$sample\t/g" | tr -d "[]()'" | sed -r "s/, /\t/g"; done > ../tables/freyja_outs.tsv
 
