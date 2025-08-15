#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# K8s layout
: "${DATA_ROOT:=/cmbg/Exom-ann/rodiny}"
: "${RUNS_ROOT:=/cmbg/Exom/data}"
: "${REF_ROOT:=/mnt/references}"
: "${PROJ_ROOT:=/cmbg/Exom-ann/src/project/kamilareblova/exomy}"
: "${METARNN:=/cmbg/Exom-ann/src/rodiny-scripts/MetaRNN}"

# VSTUPY RODINY
source ./input-rodina 

# jméno rodiny – když není dané, vezmi název adresáře
: "${jmeno_rodiny:=$(basename "$PWD")}"
cesta3="$DATA_ROOT/$jmeno_rodiny"
mkdir -p "$cesta3/bams"

# Nástroje z image
GATK="${GATK:-gatk}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
TABIX="${TABIX:-tabix}"
BGZIP="${BGZIP:-bgzip}"
SAMTOOLS="${SAMTOOLS:-samtools}"
VEP="${VEP:-vep}"
GUNZIP="${GUNZIP:-gunzip}"
PYTHON="${PYTHON:-python}"
BIOPET="${BIOPET:-java -Xmx8g -jar /cmbg/Exom-ann/src/biopet.jar}"

# rychlá kontrola dostupnosti (nevypíná běh, jen varuje)
for t in "$GATK" "$BCFTOOLS" "$TABIX" "$BGZIP" "$SAMTOOLS" "$VEP" "$GUNZIP" "$PYTHON"; do
  command -v "${t%% *}" >/dev/null || echo "[WARN] Nenalezeno v PATH: $t" >&2
done
command -v java >/dev/null || echo "[WARN] Java není v PATH – Biopet nepoběží." >&2
[[ -r "/cmbg/Exom-ann/src/biopet.jar" ]] || echo "[WARN] Chybí /cmbg/Exom-ann/src/biopet.jar" >&2


# Reference/VEP cesty (pod /mnt)
FASTA="${FASTA:-$REF_ROOT/Homo_sapiens/GRCh38/full_only/GRCh38-p10.fa}"
VEP_DIR="${VEP_DIR:-$REF_ROOT/Homo_sapiens/GATK/GRCh38/Annotation/vep}"
VEP_PLUGINS="${VEP_PLUGINS:-$VEP_DIR/Plugins}"

# --- BAMy z rodinné složky ----------------------------------------------------
# Jak vypadají jména BAMů
: "${BAM_SUFFIX:=.mdup.sorted.bam}"

BAMS_DIR="$cesta3/bams"

# najdi první BAM, který začíná ID člena (podle input-rodina)
pick_bam() { ls -1 "$BAMS_DIR/${1}"*.bam 2>/dev/null | head -n1; }

# načti podle $clen1/$clen2/$clen3
bam1="$(pick_bam "$clen1" || true)"
bam2="$(pick_bam "$clen2" || true)"
bam3="$(pick_bam "$clen3" || true)"

# kontrola, že existují
[[ -n "${bam1:-}" ]] || { echo "Chybí BAM pro $clen1 v $BAMS_DIR"; exit 2; }
[[ -n "${bam2:-}" ]] || { echo "Chybí BAM pro $clen2 v $BAMS_DIR"; exit 2; }
[[ -n "${bam3:-}" ]] || { echo "Chybí BAM pro $clen3 v $BAMS_DIR"; exit 2; }

# ------------------------------------------------------------------------------

### CONDA
# Conda z miniconda3
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Kompatibilita se starým zápisem "activate"/"deactivate"
activate()   { conda activate "$@"; }
deactivate() { conda deactivate; }


#source activate kamila
#source input-rodina
###################################

#ref="/mnt/hdd2/reference/GRCh37/full_only/GRCh37.fa"
#ref="/mnt/hdd2/reference/kristina/GRCh38/full_only/GRCh38-p10.fa"
ref="$REF_ROOT/Homo_sapiens/GRCh38/full_only/GRCh38-p10.fa"
#cesta4="/mnt/hdd2/reference"
cesta4="$REF_ROOT"

#GATK="$HOME/SOFT/gatk/gatk"
#VEP="$HOME/SOFT/vep/ensembl-vep/vep"


############################ zde napsat cestu ke vzorkum a bed filu ################
#cesta3="$PWD"
bedfile="$PROJ_ROOT/beds/HyperExomeV2_primary_targets+-50-merged.bed"

####################################################################################

# zde uvest nazvy vzorku rodiny #########################

###### Haplotypecaller
    "$GATK" --java-options "-Xmx4g" HaplotypeCaller -R $ref \
                      -I "$bam1" -I "$bam2" -I "$bam3" \
            -L $bedfile \
             --dont-use-soft-clipped-bases true \
            -A StrandBiasBySample \
             --minimum-mapping-quality 0 \
             --mapping-quality-threshold-for-genotyping 0 \
             --enable-dynamic-read-disqualification-for-genotyping true \
             --flow-filter-alleles-qual-threshold 0 \
            -O $cesta3/$jmeno_rodiny.vcf
   

"$BCFTOOLS" norm -m -both -f $ref -o $cesta3/$jmeno_rodiny.norm.vcf $cesta3/$jmeno_rodiny.vcf

######### ACGT
"$GATK" --java-options "-Xmx4g"  VariantAnnotator \
      -V $cesta3/$jmeno_rodiny.norm.vcf \
      -O $cesta3/$jmeno_rodiny.norm.acgt.vcf.gz \
      --resource:ACGT "$REF_ROOT/Homo_sapiens/ACGT/ACGT.ACCounts.noSamples.vcf.gz" \
      --expression ACGT.AF \
      --expression ACGT.AC \
      --expression ACGT.AC_Hom \
      --expression ACGT.AC_Het \
      --expression ACGT.AC_Hemi


"$GUNZIP" $cesta3/$jmeno_rodiny.norm.acgt.vcf.gz

######## VAF
"$BCFTOOLS" norm -m -both -f $ref -o $cesta3/$jmeno_rodiny.norm.vcf $cesta3/$jmeno_rodiny.norm.acgt.vcf

"$BCFTOOLS" view -O b -o $cesta3/$jmeno_rodiny.pom.bcf $cesta3/$jmeno_rodiny.norm.vcf
"$BCFTOOLS" +fill-tags $cesta3/$jmeno_rodiny.pom.bcf -Ob -o $cesta3/$jmeno_rodiny.pom2.bcf -- -t FORMAT/VAF
"$BCFTOOLS" convert -O v -o $cesta3/$jmeno_rodiny.norm.vaf.vcf $cesta3/$jmeno_rodiny.pom2.bcf

"$BCFTOOLS" annotate -a $PROJ_ROOT/beds/AR.fin.bed.gz \
    -h $PROJ_ROOT/beds/header-AR.hdr \
    -c CHROM,FROM,TO,dedicnostAR \
    -l dedicnostAR:append \
    -m -xx $cesta3/$jmeno_rodiny.norm.vaf.vcf > $cesta3/$jmeno_rodiny.pom3

"$BCFTOOLS" annotate -a $PROJ_ROOT/beds/AD.fin.bed.gz \
    -h $PROJ_ROOT/beds/header-AD.hdr \
    -c CHROM,FROM,TO,dedicnostAD  \
    -l dedicnostAD:append \
    -m -aa $cesta3/$jmeno_rodiny.pom3 > $cesta3/$jmeno_rodiny.pom4

"$BCFTOOLS" annotate -a $PROJ_ROOT/beds/Xlinked.fin.bed.gz \
    -h $PROJ_ROOT/beds/header-Xlinked.hdr \
    -c CHROM,FROM,TO,dedicnostXlinked \
    -l dedicnostXlinked:append \
    -m -bb $cesta3/$jmeno_rodiny.pom4 > $cesta3/$jmeno_rodiny.pom5

"$BCFTOOLS" annotate -a $PROJ_ROOT/beds/Ylinked.fin.bed.gz  \
    -h $PROJ_ROOT/beds/header-Ylinked.hdr \
    -c CHROM,FROM,TO,dedicnostYlinked \
    -l dedicnostYlinked:append \
    -m -cc $cesta3/$jmeno_rodiny.pom5 > $cesta3/$jmeno_rodiny.pom6

"$BCFTOOLS" annotate -a $PROJ_ROOT/beds/omim-phenotyp-new-v3.bed.gz \
    -h $PROJ_ROOT/beds/header-fenotyp.hdr \
    -c CHROM,FROM,TO,fenotyp \
    -l fenotyp:append \
    -m -yy $cesta3/$jmeno_rodiny.pom6 | sed  's/xx/dedicnostAR=0/' | grep \
    -v 'INFO=<ID=dedicnostAR=0' | sed  's/yy/fenotyp=0/' | grep \
    -v 'INFO=<ID=fenotyp=0' | sed  's/aa/dedicnostAD=0/' | grep \
    -v 'INFO=<ID=dedicnostAD=0' | sed  's/bb/dedicnostXlinked=0/' | grep \
    -v 'INFO=<ID=dedicnostXlinked=0' | sed  's/cc/dedicnostYlinked=0/' | grep \
    -v 'INFO=<ID=dedicnostYlinked=0' > $cesta3/$jmeno_rodiny.norm.vaf.omim.vcf

"$BGZIP" $cesta3/$jmeno_rodiny.norm.vaf.omim.vcf
"$TABIX" $cesta3/$jmeno_rodiny.norm.vaf.omim.vcf.gz

rm $cesta3/$jmeno_rodiny.pom*

#source activate MetaRNN
conda activate MetaRNN

"$BGZIP" $cesta3/$jmeno_rodiny.norm.vcf
"$PYTHON" $METARNN/MetaRNN.py hg38 $cesta3/$jmeno_rodiny.norm.vcf.gz
conda deactivate

rm $cesta3/*temp*

cat $cesta3/$jmeno_rodiny.norm.vcf.gz.indel.annotated $cesta3/$jmeno_rodiny.norm.vcf.gz.nsSNV.annotated > $cesta3/$jmeno_rodiny.MetaRNN

grep -v "#" $cesta3/$jmeno_rodiny.MetaRNN  > a
sed -i 's/ENST/Varsome=ENST/' a
        sed -i 's/;/ /g' a
        sed -i 's/ /*/g' a
        awk -F "\t" '{print $1, $2, ". " $3, $4,  "340", "PASS", $5, $6}' a > aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /*/g' aa
        sort -k1,1 -k2,2n aa > aaa
        cat $METARNN/hlavicka  aaa > $cesta3/$jmeno_rodiny.MetaRNN.pom.vcf

"$BGZIP" -@ 8 -f $cesta3/$jmeno_rodiny.MetaRNN.pom.vcf
"$TABIX" -p vcf $cesta3/$jmeno_rodiny.MetaRNN.pom.vcf.gz

"$GATK" --java-options "-Xmx4g"  VariantAnnotator \
     -V $cesta3/$jmeno_rodiny.norm.vaf.omim.vcf.gz \
     -O $cesta3/$jmeno_rodiny.norm.metarnn.vcf \
     --resource:MetaRNN $cesta3/$jmeno_rodiny.MetaRNN.pom.vcf.gz \
     --expression MetaRNN.Varsome

"$BGZIP" $cesta3/$jmeno_rodiny.norm.metarnn.vcf
"$TABIX" $cesta3/$jmeno_rodiny.norm.metarnn.vcf.gz

rm a aa aaa

##########ANNOVAR
$HOME/SOFT/annovar/table_annovar.pl -vcfinput "$cesta3/$jmeno_rodiny.norm.metarnn.vcf.gz" \
        "$REF_ROOT/Homo_sapiens/annovar/humandb_hg38/" \
        -buildver hg38 \
        -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20240917,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f \
        -nastring . -otherinfo \
        -polish \
        -xreffile $cesta4/gene_fullxref.txt \
        -arg '-splicing 5 -exonicsplicing',,,,,,,,,,, \
        --remove

"$BGZIP" -@ 8 -f $cesta3/$jmeno_rodiny.norm.metarnn.vcf.gz.hg38_multianno.vcf
"$TABIX" -p vcf $cesta3/$jmeno_rodiny.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz

"$GATK" --java-options "-Xmx4g" VariantsToTable -R $ref \
     -V $cesta3/$jmeno_rodiny.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz \
     -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF SB -GF VAF \
     -F dedicnostAR -F dedicnostAD -F dedicnostXlinked -F dedicnostYlinked -F fenotyp \
     -F ACGT.AF -F ACGT.AC -F ACGT.AC_Hom -F ACGT.AC_Het -F ACGT.AC_Hemi \
     -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer \
     -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all \
     -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF \
     -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F MetaRNN.Varsome -F SIFT_pred \
     -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred \
     -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer \
     -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer \
     --output $cesta3/$jmeno_rodiny.txt


############# VEP
#~/ensembl-vep/vep -i $cesta3/$jmeno_rodiny.norm.vcf.gz --cache --cache_version 114 --dir_cache /mnt/hdd2/reference/kristina/TP53-nextflow-references/Homo_sapiens/GATK/GRCh38/Annotation/vep --fasta /mnt/hdd2/reference/kristina/GRCh38/full_only/GRCh38-p10.fa --merged --offline --vcf --hgvs -o $cesta3/$jmeno_rodiny.vep.vcf --dir_plugins /mnt/hdd2/reference/kristina/TP53-nextflow-references/Homo_sapiens/GATK/GRCh38/Annotation/vep/Plugins --force_overwrite --no_stats --plugin AlphaMissense,file=/mnt/hdd2/reference/kristina/TP53-nextflow-references/Homo_sapiens/GATK/GRCh38/Annotation/vep/Plugins/AlphaMissense/AlphaMissense_hg38.tsv.gz
"$VEP" -i $cesta3/$jmeno_rodiny.norm.vcf.gz \
    --cache --cache_version 114 \
    --dir_cache "$REF_ROOT/Homo_sapiens/GATK/GRCh38/Annotation/vep" \
    --fasta "$REF_ROOT/Homo_sapiens/GRCh38/full_only/GRCh38-p10.fa" \
    --merged --offline --vcf --hgvs \
    -o $cesta3/$jmeno_rodiny.vep.vcf \
    --dir_plugins "$REF_ROOT/Homo_sapiens/GATK/GRCh38/Annotation/vep/Plugins" \
    --force_overwrite --no_stats \
    --plugin AlphaMissense,file="$REF_ROOT/Homo_sapiens/GATK/GRCh38/Annotation/vep/Plugins/AlphaMissense/AlphaMissense_hg38.tsv.gz"


#source activate biopet
#if command -v "$BIOPET" >/dev/null; then
#  "$BIOPET" tool VepNormalizer -I $cesta3/$jmeno_rodiny.vep.vcf -O $cesta3/$jmeno_rodiny.vepalpfa.anot.vcf -m standard
#else
#  echo "[WARN] Biopet není dostupný → krok VepNormalizer přeskočen." >&2
#fi

"$BIOPET" tool VepNormalizer \
  -I "$cesta3/$jmeno_rodiny.vep.vcf" \
  -O "$cesta3/$jmeno_rodiny.vepalpfa.anot.vcf" \
  -m standard

#source activate py36
conda activate py36
# python /home/kristina/TEST-VEP/vcf-simplify.py SimplifyVCF -toType table -inVCF $cesta3/$jmeno_rodiny.vepalpfa.anot.vcf -out  $cesta3/$jmeno_rodiny.vepalpfa.anot.txt
"$PYTHON" $PROJ_ROOT/scripts/vcf-simplify.py SimplifyVCF \
    -toType table \
    -inVCF $cesta3/$jmeno_rodiny.vepalpfa.anot.vcf \
    -out  $cesta3/$jmeno_rodiny.vepalpfa.anot.txt
conda deactivate

awk '{print $1, $2, $4, $5, $58, $59}' $cesta3/$jmeno_rodiny.vepalpfa.anot.txt > vyber
sed -i 's/ /\t/g' vyber
paste $cesta3/$jmeno_rodiny.txt vyber > $cesta3/$jmeno_rodiny.merged.txt



