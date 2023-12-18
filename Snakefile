[SMPS] = glob_wildcards("data/fastq/{sample}.fastq.gz")
SAMPLES = set(list(SMPS))

rule all_fastP:
   input:expand("data/trimming/{sample}.fastq.gz", sample=SMPS) 
   
rule fastP:
   input:"data/fastq/{sample}.fastq.gz",
   output:"data/trimming/{sample}.fastq.gz",
   conda:
     "env.yaml",
   shell:
     "fastp --trim_poly_x --qualified_quality_phred 20 --length_required 20 -o {output} -i {input}"

rule all_bwa:
   input:expand("data/temp/{sample}_aln.sai", sample=SMPS),
  
rule bwa_align_temp:
   input:
     "data/trimming/{sample}.fastq.gz",
   output:
     "data/temp/{sample}_aln.sai",
   conda:
     "/home/tran.986/Snakemake/env.yaml",
   shell:
     "bwa aln /fs/project/bradley.720/db/nonprokaryotic_genomes/contamination_screen/meta_contamination.fa.gz {input} > {output}"

rule all_samse:
   input:expand("data/qc/univec_core/{sample}.sam", sample=SMPS),

rule bwa_samse:
   input:
     "data/temp/{sample}_aln.sai",
     "data/trimming/{sample}.fastq.gz",
   output:
     "data/qc/univec_core/{sample}.sam",
   conda:
     "/home/tran.986/Snakemake/env.yaml",
   shell:
     "bwa samse /fs/project/bradley.720/db/nonprokaryotic_genomes/contamination_screen/meta_contamination.fa.gz {input} > {output}"
     
rule all_samtools:
    input:expand("data/filtering/contamination/{sample}.bam",sample=SMPS),

rule samtools:
    input:"data/qc/univec_core/{sample}.sam",
    output:"data/filtering/contamination/{sample}.bam",
    conda:
      "/home/tran.986/Snakemake/env.yaml",
    shell:
      "samtools view -hbS -F4 {input} > {output}"

rule all_picard:
    input:expand("data/filtering/fastq/{sample}.fastq",sample=SMPS),

rule picard:
    input:"data/filtering/contamination/{sample}.bam",
    output:"data/filtering/fastq/{sample}.fastq"
    conda:
      "env.yaml",
    shell:
      "picard SamToFastq --INPUT {input} --FASTQ {output} --VALIDATION_STRINGENCY LENIENT"
      
rule all_gunzip:
    input:expand("data/filtering/fastq/{sample}.fastq.gz",sample=SMPS),
    
rule gunzip:
    input:"data/filtering/fastq/{sample}.fastq",
    output:"data/filtering/fastq/{sample}.fastq.gz"
    conda:
      "env.yaml",
    wrapper:
      "v3.0.4/bio/bgzip"
      
rule all_megahit:
    input:
      expand("data/assembly/{sample}", sample=SMPS),
     
rule megahit:
    input:"data/filtering/fastq/{sample}.fastq.gz",
    output:
      directory("data/assembly/{sample}"),
    conda:
      "env.yaml",
    shell:
      "megahit -r {input} -o {output}"

rule all_maxbin:
     input:expand("data/filtering/fastq/{sample}.fastq.gz", sample=SMPS)

rule maxbin:
     input:
      in1="data/assembly/{sample}/final.contigs.fa",
      in2="data/filtering/fastq/{sample}.fastq.gz",
     output: 
      "data/binning/{sample}/",
     conda:
      "env.yaml"
     shell:
      "run_MaxBin.pl -contig {input.in1} -out {output} -reads {input.in2}"

rule all_checkm:
    input:[f"data/binning/{sample}/" for sample in SMPS]

rule checkm:
    input:"data/binning/{sample}/",
    output:"data/report/checkm"
    conda:
       "env.yaml",
    shell:
       "checkm taxonomy_wf domain Bacteria -x fasta {input} {output}"

rule checkm_qa:
    input:
       in1="data/reports/checkm/Bacteria.ms",
       in2="data/reports/checkm/",
    output:"data/reports/checkm/quality_JP4D.tsv",
    shell:"checkm qa {input.in1} {input.in2} --file {output} --tab_table -o 2"

#######

all_krona_cut:
    input:expand("data/taxonomy/kraken/{sample}.krona", sample=SMPS)

rule kraken:
    input:
      "data/assembly/{sample}/final.contigs.fa",
    output:
      out1="data/taxonomy/kraken/{sample}.100.kraken2",
      out2="data/taxonomy/kraken/{sample}.100.kreport",
    resources:
       time="2:00:00",
       nodes=2,
       mem_mb=60000,
    conda:"env.yaml",
    threads: 8
    shell:
       "kraken --db /fs/project/bradley.720/db/kraken_dbs/ --threads {threads} --output {output.out1} --report {output.out2}  --input {input}"
    
rule kraken_krona:
    input:"data/taxonomy/kraken/{sample}.100.kreport",
    output:"data/taxonomy/kraken/{sample}.100.krona",
    resources:
       time="2:00:00",
       nodes=2,
       mem_mb=60000,
    conda:"env.yaml",
    shell:"kreport2krona.py -r {input} -o {output}"

rule cut:
    input:
      in1="data/taxonomy/kraken/{sample}.100.kraken",
      in2="data/taxonomy/kraken/{sample}.100.krona"
    shell:"cut -f2,3 {input.in1} > {input.in2}"

rule all_ktimport:
    input:expand("data/reports/krona/{sample}.100.html", sample=SMPS)
    
rule ktimport:
    input:"data/taxonomy/kraken/{sample}.100.krona",
    output:"data/reports/krona/{sample}.100.html"
    resources:
       time="2:00:00",
       nodes=2,
       mem_mb=60000,
    conda:"env.yaml",
    shell:
       "ktImportTaxonomy {input} -o {output}

rule all_bracken:
    input:expand("data/taxonomy/bracken/{sample}.bracken", sample=SMPS)

rule bracken:
    input:"data/taxonomy/kraken/{sample}.100.kreport",
    output:
       out1="data/taxonomy/bracken/{sample}.bracken",
       out2="data/taxonomy/bracken/{sample}.reports"
    shell:
       "bracken -d /fs/project/bradley.720/db/kraken_dbs/ -i {input} -o {output.out1} -r 150 -l S -w {output.out2}


