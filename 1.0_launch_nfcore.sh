#!/bin/bash
#SBATCH --job-name=test_nextflow_all
#SBATCH --output=test_nextflow_all.out
#SBATCH --error=test_nextflow_all.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --mail-user=s.einspanier@phytomed.uni-kiel.de

export http_proxy=http://10.0.7.235:3128
export https_proxy=http://10.0.7.235:3128
export ftp_proxy=http://10.0.7.235:3128

module load gcc12-env/12.3.0
module load miniconda3/23.5.2
module load singularity
#conda activate nextflow

cd /gxfs_work/cau/suaph281/RNAseq/spen_lag_timeseries/data/0_raw_reads/




#3′ digital gene expression assays
#Some bulk RNA-seq library preparation protocols capture only a 3’ tag from each transcript, e.g. 3’Pool-seq, DRUG-seq, 
# BRB-seq or Lexogen’s commercial QuantSeq 3’ mRNA-seq FWD protocol. The following parameters have been validated for 
# QuantSeq 3' mRNA-seq FWD data, and provide useful starting points for other 3’ RNA-seq protocols:
#Custom STAR parameters

#Lexogen provides an example analysis workflow on their website, which includes the ENCODE standard 
#options for the STAR aligner. In addition, Lexogen also decreases the tolerance for mismatches and 
#clips poly(A) tails. To apply these settings, add the following parameters when running the pipeline:


        #--gtf /gxfs_work/cau/suaph281/RNAseq/references/spen/spen_GeneExt.gtf \
~/programs/nextflow/nextflow run nf-core/rnaseq \
        --input  /gxfs_home/cau/suaph281/spen_lag_phase_rnaseq/data/sample_sheet.csv \
        --outdir   /gxfs_work/cau/suaph281/RNAseq/spen_lag_timeseries/data/2024_11_12_nfcore_out/ \
        --aligner star_salmon \
        --save_unaligned \
        --fasta  /gxfs_work/cau/suaph281/RNAseq/references/spen/Spenn_top12.fasta   \
        --gtf /gxfs_work/cau/suaph281/RNAseq/references/spen/spen_GeneExt_orphan_corrected.gtf \
        --skip_biotype_qc \
        --skip_gtf_transcript_filter \
        -profile singularity \
        --extra_trimgalore_args="-a A{10}" \
        --extra_salmon_quant_args="--noLengthCorrection" \
        --extra_star_align_args "--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA" \
        --with_umi \
        --umitools_extract_method regex \
        --umitools_bc_pattern "^(?P<umi_1>.{6})(?P<discard_1>.{4}).*" 

jobinfo