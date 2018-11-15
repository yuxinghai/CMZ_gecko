shell.prefix("PATH=/home/yuxinghai/.conda/envs/py35/bin:$PATH;")
import os
import pandas as pd
base = "/data1/zhoulab/yuxinghai/CMZ_lab/GuoD"
proj = "1st"
proj_dir = os.path.join(base, proj)
raw_dir = os.path.join(proj_dir, "data", "01.cleandata")
work_dir = os.path.join(proj_dir, "results")
script_dir = os.path.join(proj_dir, "scripts")
sampName = os.path.join(base, "anno", "GuoD_CRISPRsamples_info.csv")
adapter_5 = "GTGGAAAGGACGAAACACCG"
adapter_3 = "GTTTTAGAG"
fuzzy_adapter3 = "TTTAGAG"
pool_dir = "/data2/zhoulab/yuxinghai/genome/CRISPRscreen"
idx_dir = pool_dir + "/../bt.idx/CRISPRscreen"
sample_info = pd.read_csv(sampName)
sample = sample_info.samples.unique()
cond = sample_info[sample_info.contrast != "None"].condition.unique()


rule all:
    input:
        expand("%s/fastqc/{sample}.read1_Clean_fastqc.html" % work_dir, sample = sample),
        expand("%s/filter_fq/{sample}_endBar.fq.gz" % work_dir, sample = sample),
        expand("%s/exact_sgRNA/sg_{sample}.fq" % work_dir, sample = sample),
        expand("%s/bowtie.A_idx/crisprpool_A.1.ebwt" % idx_dir, sample = sample),
        expand("%s/map/{sample}_unB.fq" % work_dir, sample = sample),
        expand("%s/count/{sample}.screen_A.csv" % work_dir, sample = sample),
        expand("%s/run_qc/{condition}_cdf.pdf" % work_dir, condition = cond),
        expand("%s/run_mageck/{condition}/{condition}.A.tsv" % work_dir, condition = cond),


rule fastqc:
    input:
        R1 = raw_dir + "/{sample}.read1_Clean.fastq.gz",
        R2 = raw_dir + "/{sample}.read2_Clean.fastq.gz",
    output:
        "%s/fastqc/{sample}.read1_Clean_fastqc.html" % work_dir,
        "%s/fastqc/{sample}.read2_Clean_fastqc.html" % work_dir,
    threads: 2
    shell:
        '''
	    fastqc -f fastq -t {threads} --noextract -o {work_dir}/fastqc {input.R1} {input.R2}
	    '''

rule filter_fq:
    input:
        R1 = raw_dir + "/{sample}.read1_Clean.fastq.gz",
        R2 = raw_dir + "/{sample}.read2_Clean.fastq.gz",
    output:
        filter_fq = temp("%s/filter_fq/{sample}_filter.fq.gz" % work_dir),
        endBar_fq = "%s/filter_fq/{sample}_endBar.fq.gz" % work_dir,
        stat = "%s/filter_fq/{sample}.stat" % work_dir,
    shell:
        '''
        zcat {input.R1} {input.R2} | paste -d "\!" - - - - | grep "{fuzzy_adapter3}" | \
        tr "\!" "\\n" | gzip -c > {output.filter_fq}
        python {script_dir}/filter_fq.py {output.filter_fq} {output.endBar_fq}
        wc -l {work_dir}/filter_fq/{wildcards.sample}* > {output.stat}
        '''

rule exact_sgRNA:
    input:
        endBar_fq = "%s/filter_fq/{sample}_endBar.fq.gz" % work_dir,
    output:
        sg_fq = "%s/exact_sgRNA/sg_{sample}.fq" % work_dir,
    threads: 5
    log:
        "%s/exact_sgRNA/{sample}.log" % work_dir,
    shell:
        '''
        cutadapt -g "{adapter_5};anywhere" -a "{adapter_3};anywhere" \
        -e 0.1 -j {threads} -m 17 -M 23 -n 2 -o {output.sg_fq} {input.endBar_fq} 1>{log} 2>&1
        '''

rule build_index:
    input:
        csv_a = pool_dir + "/human_geckov2_library_a.csv",
        csv_b = pool_dir + "/human_geckov2_library_b.csv",
    output:
        fa_a = pool_dir + "/crisprpool_A.fa",
        fa_b = pool_dir + "/crisprpool_B.fa",
        ebwt_A = "%s/bowtie.A_idx/crisprpool_A.1.ebwt" % idx_dir,
        ebwt_B = "%s/bowtie.B_idx/crisprpool_B.1.ebwt" % idx_dir,
    shell:
        """
        python {script_dir}/mkpool.py {pool_dir} {input.csv_a} {input.csv_b}
        bowtie-build {output.fa_a} {idx_dir}/bowtie.A_idx/crisprpool_A
        bowtie-build {output.fa_b} {idx_dir}/bowtie.B_idx/crisprpool_B
        """

rule map:
    input:
        sg_fq = "%s/exact_sgRNA/sg_{sample}.fq" % work_dir,
        ebwt_A = "%s/bowtie.A_idx/crisprpool_A.1.ebwt" % idx_dir,
        ebwt_B = "%s/bowtie.B_idx/crisprpool_B.1.ebwt" % idx_dir,
    output:
        A_sam = temp("%s/map/{sample}.A_sam" % work_dir),
        B_sam = temp("%s/map/{sample}.B_sam" % work_dir),
        A_bam = "%s/map/{sample}.A_bam" % work_dir,
        B_bam = "%s/map/{sample}.B_bam" % work_dir,
        A_un = "%s/map/{sample}_unA.fq" % work_dir,
        B_un = "%s/map/{sample}_unB.fq" % work_dir,
    threads: 12
    log:
        A = "%s/map/{sample}.Alog" % work_dir,
        B = "%s/map/{sample}.Blog" % work_dir,
    shell:
        """
        bowtie -p {threads} --best --strata --norc -m 1 -l 17 -n 3 \
        --un {output.A_un} -S {idx_dir}/bowtie.A_idx/crisprpool_A  {input.sg_fq} {output.A_sam} 1>{log.A} 2>&1
        samtools view -@ {threads} -bS {output.A_sam} > {output.A_bam}
        bowtie -p {threads} --best --strata --norc -m 1 -l 17 -n 3 \
        --un {output.B_un} -S {idx_dir}/bowtie.B_idx/crisprpool_B {input.sg_fq} {output.B_sam} 1>{log.B} 2>&1
        samtools view -@ {threads} -bS {output.B_sam} > {output.B_bam}
        """

rule count:
    input:
        A_bam = "%s/map/{sample}.A_bam" % work_dir,
        csv_a = pool_dir + "/human_geckov2_library_a.csv",
        csv_b = pool_dir + "/human_geckov2_library_b.csv",
    output:
        csv_A = "%s/count/{sample}.screen_A.csv" % work_dir,
    log:
        "%s/count/{sample}.log" % work_dir,
    shell:
        """
        python {script_dir}/shcount.py {input.csv_a} {input.A_bam} -o {output.csv_A} 1>{log} 2>&1
        """

rule run_qc:
    input:
        contrast = lambda wildcards : [work_dir + "/count/" + str(sap) +  ".screen_A.csv" 
            for sap in sample_info[sample_info.condition == wildcards.condition].contrast.unique()],
        treat = lambda wildcards :[work_dir + "/count/" + str(sap) +  ".screen_A.csv"
            for sap in sample_info[sample_info.condition == wildcards.condition].samples.unique()],
        Lib_A = "%s/count/Lib-A.screen_A.csv" % work_dir,
    output:
        cdf = "%s/run_qc/{condition}_cdf.pdf" % work_dir,
        skew = "%s/run_qc/{condition}_skew.tsv" % work_dir,
    run:
        str_csv = ":".join(input)
        qc_script = script_dir + "/qc.R"
        command = " ".join(["Rscript", qc_script, work_dir, str_csv, output.cdf, output.skew])
        print(command)
        shell(command)

rule run_mageck:
    input:
        contrast = lambda wildcards : [work_dir + "/count/" + str(sap) +  ".screen_A.csv" 
            for sap in sample_info[sample_info.condition == wildcards.condition].contrast.unique()],
        condition = lambda wildcards :[work_dir + "/count/" + str(sap) +  ".screen_A.csv"
            for sap in sample_info[sample_info.condition == wildcards.condition].samples.unique()],
        sampName = sampName           
    output:
        count = "%s/run_mageck/{condition}/{condition}.A.tsv" % work_dir,
        sgrna = "%s/run_mageck/{condition}/{condition}.sgrna_summary.txt" % work_dir,
        gene = "%s/run_mageck/{condition}/{condition}.gene_summary.txt" % work_dir,
    run:
        contrast = sample_info[sample_info.condition == wildcards.condition].contrast.unique()[0]
        str_con = ":".join([wildcards.condition, contrast])
        R_script = script_dir + "/pre_Mageck.R"
        Rcommand = " ".join(["Rscript", R_script, work_dir, str_con, sampName, output.count])
        prefix = os.path.join(work_dir, "run_mageck", wildcards.condition, wildcards.condition)
        control = ",".join(sample_info[sample_info.condition == wildcards.condition].contrast.unique())
        treat = ",".join(sample_info[sample_info.condition == wildcards.condition].samples.unique())
        Mcommand = " ".join(["mageck test -k", output.count, "-t", treat, "-c", control,\
            "--norm-method total --normcounts-to-file --sort-criteria pos --adjust-method fdr -n",\
            prefix, "--pdf-report"])
        print(Rcommand)
        shell(Rcommand)
        print(Mcommand)
        shell(Mcommand)