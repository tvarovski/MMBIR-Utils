configfile: "/Users/temagbetere/workspace/nfsscratchTools_shared/mmbir_pipeline_config_mouse.yaml"


rule all:
    input:
        "{bamin}.counts.txt",
        submit_arr=expand("{base_path}/clusterArrayT.{{bamin}}.sh", base_path=config["outputs_home"])

rule removeDups:
    params:
        script=config["paths"]["scripts"]["markdups"]
    threads: config["threads"]
    input:
        "{bamin}.bam"
    output:
        "{bamin}.deduped.bam",
        "{bamin}.dedup_metrics.txt"
    shell:
        "{params.script} {input}"

rule getUnaligned:
    params:
        script=config["paths"]["scripts"]["getUnaligned"]
    threads: config["threads"]
    input:
        "{bamin}.deduped.bam"
    output:
        "{bamin}.UM_bam",
        "{bamin}.UM_fq"
    shell:
        "{params.script} {input} {output}"

rule getIndels:
    params:
        script=config["paths"]["scripts"]["getInDels"]
    threads: config["threads"]
    input:
        "{bamin}.deduped.bam"
    output:
        "{bamin}.idsc_bam",
        "{bamin}.idsc_fq"
    shell:
        "{params.script} {input} {output}"

rule realignUnaligned:
    params:
        script=config["paths"]["scripts"]["realign_UM"]
    threads: config["threads"]
    input:
        "{bamin}.UM_fq"
    output:
        "{bamin}.real_bam"
    shell:
        "{params.script} {input} {output}"

rule getRealignedUM:
    params:
        script=config["paths"]["scripts"]["getUnaligned"]
    threads: config["threads"]
    input:
        "{bamin}.real_bam"
    output:
        "{bamin}.real_UM_fq",
        "{bamin}.real_UM_bam"
    shell:
        "{params.script} {input} {output}"

rule getRealignedIndels:
    params:
        script=config["paths"]["scripts"]["getInDels"]
    threads: config["threads"]
    input:
        "{bamin}.real_bam"
    output:
        "{bamin}.real_idsc_bam",
        "{bamin}.real_idsc_fq"
    shell:
        "{params.script} {input} {output}"

rule create_UM_idsc_fq:
    threads: config["threads"]
    input:
        "{bamin}.UM_fq",
        "{bamin}.idsc_fq",
        "{bamin}.real_UM_fq",
        "{bamin}.real_idsc_fq"
    output:
        "{bamin}.fq"
    shell:
        "cat {input} > {output}"


rule filter_reads_trim:
    params:
        trimmomatic=config["paths"]["trimmomatic_home"],
        adapters=config["paths"]["adapters_fasta"]
    threads: config["threads"]
    input:
        "{bamin}.fq"
    output:
        "{bamin}.fq.pp"
    shell:
        "module unload jdk\nmodule load stack/legacy\nmodule load openjdk/11.0.2\njava -jar {params.trimmomatic} SE -threads {threads} -phred33 {input} {output} ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule filter_reads_prinseq:
    threads: config["threads"]
    input:
        "{bamin}.fq.pp"
    output:
        "{bamin}.PP_fq"
    shell:
        "prinseq++ -threads {threads} -fastq {input} -lc_dust=0.07 -out_good {output} -out_bad /dev/null"

rule changeFastqHeader:
    params:
        script=config["paths"]["scripts"]["changeFastqHeader"]
    input:
        "{bamin}.PP_fq"
    output:
        "{bamin}.PPh.fq"
    shell:
        "{params.script} {input} {output}"

rule count_reads:
    input:
        "{bamin}.deduped.bam",
        "{bamin}.fq"
    output:
        "{bamin}.counts.txt"
    shell:
        "echo {bamin}.deduped.bam >> counts.txt\nsamtools view -c -@ {threads} {bamin}.deduped.bam >> counts.txt\necho {bamin}.fq >> counts.txt\nsamtools view -c -@ {threads} {bamin}.fq >> counts.txt"

rule CreateDir:
    params:
        outputs_home=config["outputs_home"]
    input:
        "{bamin}.PPh.fq"
    output:
        expand("{home}/{{bamin}}/{{bamin}}.PPh.fq", home=config["outputs_home"])
    shell:
        "mkdir -p {params.outputs_home}/{wildcards.bamin}/\nmv {input} {params.outputs_home}/{wildcards.bamin}/"

rule createConfig2:
    params:
        outputs_home=config["outputs_home"],
        config2=config["paths"]["base_scripts"]["config2"]
    input:
        in_fq=expand("{outputs_home}/{{bamin}}/{{bamin}}.PPh.fq", outputs_home=config["outputs_home"])
    output:
        expand("{MMBSearch_home}/config2.{{bamin}}.txt", MMBSearch_home=config["MMBSearch_home"])
    shell:
        "sed 's/BASE/{wildcards.bamin}/g' {params.config2} > {output}"

rule mmbirAln:
    params:
        outputs_home=config["outputs_home"],
        MMBSearch_home=config["MMBSearch_home"],
        script=config["paths"]["scripts"]["MMBSearch_Aln"]
    threads: config["threads"]
    input:
        config2_file=expand("{MMBSearch_home}/config2.{{bamin}}.txt", MMBSearch_home=config["MMBSearch_home"]),
        input_fq_file=expand("{outputs_home}/{{bamin}}/{{bamin}}.PPh.fq", outputs_home=config["outputs_home"])
    output:
        expand("{outputs_home}/{{bamin}}/bwaAligned.sam", outputs_home=config["outputs_home"])
    shell:
        "cd {params.MMBSearch_home}\n{params.script} {params.MMBSearch_home}/config2.{wildcards.bamin}.txt"

rule createConfig_split:
    params:
        outputs_home=config["outputs_home"],
        config_split_base=config["paths"]["base_scripts"]["config_split"]
    input:
        expand("{outputs_home}/{{bamin}}/bwaAligned.sam", outputs_home=config["outputs_home"])
    output:
        expand("{MMBSearch_home}/config-split_{{bamin}}.txt",  MMBSearch_home=config["MMBSearch_home"])
    shell:
        "sed 's/BASE/{wildcards.bamin}/g' {params.config_split_base} > {output}"

rule splitByChr:
    params:
        MMBSearch_home=config["MMBSearch_home"],
        run_split=config["paths"]["scripts"]["run_config-split"]
    input:
        expand("{MMBSearch_home}/config-split_{{bamin}}.txt", MMBSearch_home=config["MMBSearch_home"]),
        expand("{outputs_home}/{{bamin}}/{{bamin}}.PPh.fq", outputs_home=config["outputs_home"]),
        expand("{outputs_home}/{{bamin}}/bwaAligned.sam", outputs_home=config["outputs_home"])
    output:
        expand("{outputs_home}/{{bamin}}/bwaAligned0000.sam", outputs_home=config["outputs_home"])
    shell:
        "{params.run_split} {params.MMBSearch_home}/config-split_{wildcards.bamin}.txt"

rule createChrDir:
    params:
        outputs_home=config["outputs_home"]
    input:
        expand("{outputs_home}/{{bamin}}/bwaAligned0000.sam", outputs_home=config["outputs_home"])
    output:
        expand("{outputs_home}/{{bamin}}/chr1/bwaAligned0000.sam", outputs_home=config["outputs_home"])
    shell:
        "cd {params.outputs_home}\n./makeChrDir_v2.sh {wildcards.bamin}"

rule createclusterArrs:
    params:
        script=config["paths"]["base_scripts"]["cluster_arr"]
    input:
         expand("{outputs_home}/{{bamin}}/chr1/bwaAligned0000.sam", outputs_home=config["outputs_home"])
    output:
         expand("{outputs_home}/clusterArrayT.{{bamin}}.sh", outputs_home=config["outputs_home"])
    shell:
        "sed 's/BASE/{wildcards.bamin}/g' {params.script} > {output}"

