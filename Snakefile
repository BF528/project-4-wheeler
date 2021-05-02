
SAMPLES = ["SRR3879604", "SRR3879605", "SRR3879606"]

rule all:
	input: expand("/projectnb/bf528/users/wheeler/project_4/samples/{sample}_bc_count", sample=SAMPLES)

rule count_bc_read:
        input: "/projectnb/bf528/project_4_scrnaseq/fastq/{sample}/{sample}_1_bc.fastq.gz"
        output: "/projectnb/bf528/users/wheeler/project_4/samples/{sample}_bc_count.txt"
        shell: "zcat {input} | awk 'NR%4==2' | cut -c -19 | sort | uniq -c | sort > {output}"

