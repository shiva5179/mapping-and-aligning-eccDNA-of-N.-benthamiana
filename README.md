# mapping-and-aligning-eccDNA-of-N.-benthamiana
conda create -n eccdna_env python=3.9 minimap2 samtools -y conda activate eccdna_env
chmod +x mapping-and-aligning-eccDNA-of-N.-benthamiana.py
python new\ pipline.py --input input.fastq \ --output unmapped_reads.bam \ --mito mito_ref.fasta \ --chloro chloro_ref.fasta \ --viral_adna viral_adna_ref.fasta \ --viral_bsat viral_bsat_ref.fasta \ --final_fastq final_output.fastq
