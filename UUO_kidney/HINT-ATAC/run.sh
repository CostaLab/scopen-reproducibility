#!/bin/bash
#
#### Job name
#SBATCH -J diff_footprinting
#SBATCH -e ./diff_footprinting2.txt
#SBATCH -o ./diff_footprinting2.txt
#SBATCH -t 120:00:00
#SBATCH --mem=960G --cpus-per-task=96

source ~/.bashrc
conda activate r-4.0.3

bam_loc=/data/scATA/SingleCellOpenChromatin/local/UUO_scATAC/HINT/DiffFootprintsAllCelltypes/BAM
bigwig_loc=/data/scATA/SingleCellOpenChromatin/local/UUO_scATAC/HINT/DiffFootprintsAllCelltypes/BigWig
peaks_loc=/data/scATA/SingleCellOpenChromatin/local/UUO_scATAC/HINT/DiffFootprintsAllCelltypes/Peaks
footprint_loc=/data/scATA/SingleCellOpenChromatin/local/UUO_scATAC/HINT/DiffFootprintsAllCelltypes/Footprints
motifmatching_loc=/data/scATA/SingleCellOpenChromatin/local/UUO_scATAC/HINT/DiffFootprintsAllCelltypes/MotifMatching
diff_loc=/data/scATA/SingleCellOpenChromatin/local/UUO_scATAC/HINT/DiffFootprintsAllCelltypes/Diff

mkdir -p ${bigwig_loc}
mkdir -p ${peaks_loc}
mkdir -p ${footprint_loc}
mkdir -p ${motifmatching_loc}
mkdir -p ${bam_loc}
mkdir -p ${diff_loc}

samtools merge -f --threads 96 ${bam_loc}/PT_S2.bam ../SplitBamFilesByClustering/BAM/*_1.bam
samtools merge -f --threads 96 ${bam_loc}/TAL.bam ../SplitBamFilesByClustering/BAM/*_2.bam ../SplitBamFilesByClustering/BAM/*_17.bam
samtools merge -f --threads 96 ${bam_loc}/PT_S3.bam ../SplitBamFilesByClustering/BAM/*_3.bam
samtools merge -f --threads 96 ${bam_loc}/DCT.bam ../SplitBamFilesByClustering/BAM/*_5.bam
samtools merge -f --threads 96 ${bam_loc}/PT_S1.bam ../SplitBamFilesByClustering/BAM/*_6.bam
samtools merge -f --threads 96 ${bam_loc}/CD_PC.bam ../SplitBamFilesByClustering/BAM/*_7.bam
samtools merge -f --threads 96 ${bam_loc}/EC.bam ../SplitBamFilesByClustering/BAM/*_8.bam
samtools merge -f --threads 96 ${bam_loc}/CNT.bam ../SplitBamFilesByClustering/BAM/*_9.bam
samtools merge -f --threads 96 ${bam_loc}/IC.bam ../SplitBamFilesByClustering/BAM/*_10.bam
samtools merge -f --threads 96 ${bam_loc}/Fib.bam ../SplitBamFilesByClustering/BAM/*_11.bam
samtools merge -f --threads 96 ${bam_loc}/DL_TAL.bam ../SplitBamFilesByClustering/BAM/*_12.bam
samtools merge -f --threads 96 ${bam_loc}/Mac.bam ../SplitBamFilesByClustering/BAM/*_13.bam
samtools merge -f --threads 96 ${bam_loc}/Injured_PT.bam ../SplitBamFilesByClustering/BAM/*_14.bam
samtools merge -f --threads 96 ${bam_loc}/Lymphoid.bam ../SplitBamFilesByClustering/BAM/*_15.bam
samtools merge -f --threads 96 ${bam_loc}/Pod.bam ../SplitBamFilesByClustering/BAM/*_16.bam


footprinting(){
    	filename=$1
    	samtools index $1
    	macs2 callpeak -t $filename -n ${filename%.bam} --outdir $2 -g mm --nomodel -f BAMPE -q 0.01 --keep-dup all
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_peaks.narrowPeak
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_summits.bed
    	bamCoverage -bs 24 -b $filename -o $3/${filename%.bam}.bw -p 24 --normalizeUsing CPM
	rgt-hint footprinting --organism=mm10 --atac-seq --paired-end --output-location=$4 --output-prefix=${filename%.bam} $1 $2/${filename%.bam}_peaks.narrowPeak
	rgt-motifanalysis matching --organism=mm10 --output-location=$5 --input-files $4/${filename%.bam}.bed
}

cd ${bam_loc}
for filename in *.bam;
do
	footprinting $filename ${peaks_loc} ${bigwig_loc} ${footprint_loc} ${motifmatching_loc} & 
done

wait

cd ${diff_loc}
rgt-hint differential --organism=mm10 --bc --nc 64 \
--mpbs-files=${motifmatching_loc}/CD_PC_mpbs.bed,\
${motifmatching_loc}/CNT_mpbs.bed,\
${motifmatching_loc}/DCT_mpbs.bed,\
${motifmatching_loc}/DL_TAL_mpbs.bed,\
${motifmatching_loc}/EC_mpbs.bed,\
${motifmatching_loc}/Fib_mpbs.bed,\
${motifmatching_loc}/IC_mpbs.bed,\
${motifmatching_loc}/Injured_PT_mpbs.bed,\
${motifmatching_loc}/Lymphoid_mpbs.bed,\
${motifmatching_loc}/Mac_mpbs.bed,\
${motifmatching_loc}/Pod_mpbs.bed,\
${motifmatching_loc}/PT_S1_mpbs.bed,\
${motifmatching_loc}/PT_S2_mpbs.bed,\
${motifmatching_loc}/PT_S3_mpbs.bed,\
${motifmatching_loc}/TAL_mpbs.bed \
--reads-files=${bam_loc}/CD_PC.bam,\
${bam_loc}/CNT.bam,\
${bam_loc}/DCT.bam,\
${bam_loc}/DL_TAL.bam,\
${bam_loc}/EC.bam,\
${bam_loc}/Fib.bam,\
${bam_loc}/IC.bam,\
${bam_loc}/Injured_PT.bam,\
${bam_loc}/Lymphoid.bam,\
${bam_loc}/Mac.bam,\
${bam_loc}/Pod.bam,\
${bam_loc}/PT_S1.bam,\
${bam_loc}/PT_S2.bam,\
${bam_loc}/PT_S3.bam,\
${bam_loc}/TAL.bam \
--conditions=CD_PC,CNT,DCT,DL_TAL,EC,Fib,IC,Injured_PT,Lymphoid,Mac,Pod,PT_S1,PT_S2,PT_S3,TAL \
--output-location=${diff_loc} \
--output-prefix=All





