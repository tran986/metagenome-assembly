[SMPS] = glob_wildcards("data/fastq/{sample}.fastq.gz")
SAMPLES = set(list(SMPS))

rule all_fastP:
   input:expand("data/trimming/{sample}.fastq.gz", sample=SMPS) 
   
rule fastP:
   input:"data/fastq/{sample}.fastq.gz",
   output:"data/trimming/{sample}.fastq.gz",
   resources:
     time="2:00:00",
     nodes=1,
     mem_mb=20000,
   conda:
     "env.yaml",
   shell:
     "fastp --trim_poly_x --qualified_quality_phred 20 --length_required 20 -o {output} -i {input}"

#using bwa alignment to univec_core database: 
#(fastq --> sam and )

rule all_bwa:
   input:expand("data/temp/{sample}_aln.sai", sample=SMPS),
  
rule bwa_align_temp:
   input:
     "data/trimming/{sample}.fastq.gz",
   output:
     "data/temp/{sample}_aln.sai",
   resources:
     runtime_min=300,
     nodes=1,
     mem_mb=20000,
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
   resources:
     runtime_min=300,
     nodes=2,
     mem_mb=20000,
   conda:
     "/home/tran.986/Snakemake/env.yaml",
   shell:
     "bwa samse /fs/project/bradley.720/db/nonprokaryotic_genomes/contamination_screen/meta_contamination.fa.gz {input} > {output}"
     
rule all_samtools:
    input:expand("data/filtering/contamination/{sample}.bam",sample=SMPS),

rule samtools:
    input:"data/qc/univec_core/{sample}.sam",
    output:"data/filtering/contamination/{sample}.bam",
    resources:     
      runtime_min=300,
      nodes=2,
      mem_mb=60000,
    conda:
      "/home/tran.986/Snakemake/env.yaml",
    shell:
      "samtools view -hbS -F4 {input} > {output}"

rule all_picard:
    input:expand("data/filtering/fastq/{sample}.fastq",sample=SMPS),

rule picard:
    input:"data/filtering/contamination/{sample}.bam",
    output:"data/filtering/fastq/{sample}.fastq"
    resources:
      time="2:00:00",
      nodes=2,
      mem_mb=60000,
    conda:
      "env.yaml",
    shell:
      "picard SamToFastq --INPUT {input} --FASTQ {output} --VALIDATION_STRINGENCY LENIENT"
      
rule all_gunzip:
    input:expand("data/filtering/fastq/{sample}.fastq.gz",sample=SMPS),
    
rule gunzip:
    input:"data/filtering/fastq/{sample}.fastq",
    output:"data/filtering/fastq/{sample}.fastq.gz"
    resources:
      time="2:00:00",
      nodes=2,
      mem_mb=60000,
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
    resources:
      time="2:00:00",
      nodes=2,
      mem_mb=60000,
    conda:
      "env.yaml",
    shell:
      "megahit -r {input} -o {output}"

"""
rule all_metawrap:
    input:expand("data/taxonomy/metawrap/classify_bins/{sample}/", sample=SMPS),

rule metawrap_bob:
    input:
      in1="data/assembly/{sample}/final.contigs.fa",
      in2="data/filtering/fastq/{sample}.fastq.gz",
    output:
      out1="data/reports/blobology/{sample}/",
      out2="data/binning/metawrap/binning_refinement/{sample}/",
    threads:2
    resources:
      time="2:00:00",
      nodes=2,
      mem_mb=60000,
    conda:"env_metawrap.yaml",
    shell:
      "metawrap blobology -a {input.in1} -t {threads} -o {output.out1} --bins {output.out2} {input.in2}"


rule metawrap_bin:
    input:"data/binning/metawrap/binning_refinement/{sample}/",
    output:"data/taxonomy/metawrap/classify_bins/{sample}/",
    threads:2
    resources:
      time="2:00:00",
      nodes=2,
      mem_mb=60000,
    conda:"env_metawrap.yaml",
    shell:
      "metawrap classify_bins -b {input} -o {output} -t {threads}"
"""
