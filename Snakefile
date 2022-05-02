RAWDATA_DIR = "inputs/raw"
SAMPLES=['BIS_CAN_1', 'BIS_CAN_2', 'BIS_CAN_4', 'BIS_CAN_8', 
         'BOAR_AU_13', 'BOAR_AU_4', 'BOAR_AU_6', 'BOAR_AU_8', 
         'COW_VT_2', 'COW_VT_3', 'COW_VT_6', 'COW_VT_9', 
         'dab_15', 'DAB_7', 'frida_15', 'frida_7', 
         'gp1', 'gp2', 'gp3', 'gp5', 
         'gp6_lab', 'gp7_lab', 'gp10_lab', 'gp16_lab',
         'L1_0', 'L11_0', 'L11_28', 'L16_0', 'L16_28', 'L2_0', 
         'L2_28', 'L20_0', 'L20_28', 'L4_0', 'L4_28', 'L6_0', 'L6_28', 
         'logan_14', 'logan_7', 'mardy_15', 'Mardy_7', 
         'melody_14', 'melody_7', 'nova_14', 'nova_7', 
         'PIG_VT_3', 'PIG_VT_7', 'PIG_VT_8', 'PIG_VT_9', 
         'RAB_MA_10', 'RAB_MA_3', 'RAB_MA_5', 'RAB_MA_8', 
         'RAB_PT__01', 'RAB_PT__02', 'RAB_PT__05', 'RAB_PT__12', 
         'RAT_MA_2', 'RAT_MA_3', 'RAT_MA_5', 'RAT_MA_6', 
         'RAT_NYC_1', 'RAT_NYC_2', 'RAT_NYC_5', 'RAT_NYC_6', 
         'SHE_VT_309', 'SHE_VT_402', 'SHE_VT_604', 'SHE_VT_mia',
         'SHE_WY_17_228_BHS', 'SHE_WY_17_237', 'SHE_WY_17_238', 'SHE_WY_18_053_BHS', 
         'snowman_14', 'snowman_7', 'Templeton_15', 'Templeton_7', 
         'W1_0', 'W1_28', 'W11_0', 'W11_28', 
         'W12_0', 'W12_28', 'W13_0', 'W13_28', 
         'W3_0', 'W3_28', 'W4_0', 'W4_28', 
         'wildmouse_4', 'wildmouse_5', 'wildmouse_6', 'wildmouse2', 
         'wildmouse3', 'wildmouse4', 'wildmouse7']

rule all:
    input:
        'outputs/comp/comp_nohost.csv',
        expand('outputs/map_to_megahit/{sample}.flagstat', sample = SAMPLES),
        expand("outputs/gather/{sample}.csv", sample = SAMPLES)

rule cat_samples_R1:
# concatenate lanes together into single sample
    output: 'inputs/cat/{sample}_R1.fq.gz'
    params: indir = RAWDATA_DIR
    shell:'''
    cat {params.indir}/{wildcards.sample}_S*_R1_001.fastq.gz > {output}
    '''

rule cat_samples_R2:
# concatenate lanes together into single sample
    output: 'inputs/cat/{sample}_R2.fq.gz'
    params: indir = RAWDATA_DIR
    shell:'''
    cat {params.indir}/{wildcards.sample}_S*_R2_001.fastq.gz > {output}
    '''

rule download_adapters:
    output: "inputs/adapters.fa"
    shell:'''
    # place holder to download all_PE.fa. Generated by concatenating all PE
    # illumina adapters together that ship with trimmomatic-0.39.
    '''

rule adapter_trim:
    input:
        r1 = "inputs/cat/{sample}_R1.fq.gz",
        r2 = 'inputs/cat/{sample}_R2.fq.gz',
        adapters = 'inputs/adapters.fa'
    output:
        r1 = 'outputs/trim/{sample}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{sample}_R2.trim.fq.gz',
        o1 = 'outputs/trim/{sample}_o1.trim.fq.gz',
        o2 = 'outputs/trim/{sample}_o2.trim.fq.gz'
    conda: 'trimmomatic.yml'
    shell:'''
     trimmomatic PE {input.r1} {input.r2} \
             {output.r1} {output.o1} {output.r2} {output.o2} \
             ILLUMINACLIP:{input.adapters}:2:0:15 MINLEN:25  \
             LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2
    '''
     
rule download_host:
    output:
        bison='inputs/host/bison.fna.gz',
        pig='inputs/host/pig.fna.gz',
        cow='inputs/host/cow.fna.gz',
        gp='inputs/host/gp.fna.gz',
        rabbit='inputs/host/rabbit.fna.gz',
        rat='inputs/host/rat.fna.gz',
        sheep='inputs/host/sheep.fna.gz',
        mouse='inputs/host/mouse.fna.gz',
        dog='inputs/host/dog.fna.gz',
        wolf='inputs/host/wolf.fna.gz'
    shell:'''
    wget -O {output.rat} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.fna.gz
    wget -O {output.pig} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz
    wget -O {output.sheep} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/742/125/GCF_002742125.1_Oar_rambouillet_v1.0/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.fna.gz
    wget -O {output.rabbit} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/625/GCF_000003625.3_OryCun2.0/GCF_000003625.3_OryCun2.0_genomic.fna.gz
    wget -O {output.mouse} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
    wget -O {output.gp} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/735/GCF_000151735.1_Cavpor3.0/GCF_000151735.1_Cavpor3.0_genomic.fna.gz
    wget -O {output.bison} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/754/665/GCF_000754665.1_Bison_UMD1.0/GCF_000754665.1_Bison_UMD1.0_genomic.fna.gz
    wget -O {output.cow} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz
    wget -O {output.dog} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_genomic.fna.gz
    wget -O {output.wolf} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/922/845/GCA_007922845.1_UniMelb_Wolf_Refassem_1/GCA_007922845.1_UniMelb_Wolf_Refassem_1_genomic.fna.gz
    '''
                           
rule remove_host:
    output:
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz',
        host_r1='outputs/bbduk/{sample}_R1.host.fq.gz',
        host_r2='outputs/bbduk/{sample}_R2.host.fq.gz'
    input:
        r1 = 'outputs/trim/{sample}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{sample}_R2.trim.fq.gz',
        bison='inputs/host/bison.fna.gz',
        pig='inputs/host/pig.fna.gz',
        cow='inputs/host/cow.fna.gz',
        gp='inputs/host/gp.fna.gz',
        rabbit='inputs/host/rabbit.fna.gz',
        rat='inputs/host/rat.fna.gz',
        sheep='inputs/host/sheep.fna.gz',
        mouse='inputs/host/mouse.fna.gz',
        dog='inputs/host/dog.fna.gz',
        wolf='inputs/host/wolf.fna.gz'
    run:
        wolves = ['dab_15', 'DAB_7', 'frida_15', 'frida_7', 'mardy_15', 
                  'Mardy_7', 'Templeton_15', 'Templeton_7']
        dogs = ['logan_14', 'logan_7', 'snowman_14', 'snowman_7', 
                'melody_14', 'melody_7', 'nova_14', 'nova_7']
        if "BIS_CAN" in wildcards.sample:
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.bison}")
        elif "BOAR_AU" in wildcards.sample:
            shell("echo boar {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.pig}")
        elif "COW_VT" in wildcards.sample:
            shell("echo cow {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.cow}")
        elif "gp" in wildcards.sample:
            shell("echo gp {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.gp}")
        elif "PIG_VT" in wildcards.sample:
            shell("echo pig {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.pig}")
        elif "RAB" in wildcards.sample:
            shell("echo rabbit {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.rabbit}")
        elif "RAT" in wildcards.sample:
            shell("echo rat {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.rat}")
        elif "SHE" in wildcards.sample:
            shell("echo sheep {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.sheep}")
        elif "wildmouse" in wildcards.sample:
            shell("echo wildmouse {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.mouse}")
        elif wildcards.sample in wolves:
            shell("echo wolf {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.wolf}")
        elif wildcards.sample in dogs:
            shell("echo dog {wildcards.sample}")             
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.dog}")
        else: 
            shell("echo mouse {wildcards.sample}")
            shell("bbduk.sh -Xmx64g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.host_r1} outm2={output.host_r2} k=31 ref={input.mouse}")


rule sourmash_compute:
    input:
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz'
    output: 'outputs/sigs/{sample}.sig'
    conda: 'sourmash.yml'
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance \
            --merge {wildcards.sample} -o {output} \
            {input.r1} {input.r2}
    '''

rule sourmash_compare:
    input: expand('outputs/sigs/{sample}.sig', sample = SAMPLES)
    output: 'outputs/comp/comp_nohost.csv'
    conda: 'sourmash.yml'
    shell:'''
    sourmash compare -k 31 --csv {output} {input}
    '''

rule sourmash_gather:
    input: 
        sig='outputs/sigs/{sample}.sig',
        genbank='../cosmo-kmers/inputs/databases/genbank-d2-k51.sbt.json'
    output: "outputs/gather/{sample}.csv"
    conda: 'sourmash.yml'
    shell:'''
    sourmash gather -o {output} --scaled 2000 -k 51 {input.sig} {input.genbank}
    '''
 
rule megahit:
    input:
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz'
    output: 'outputs/megahit/{sample}.contigs.fa'
    conda: 'megahit.yml'
    params: output_folder = 'outputs/megahit/'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to occur
    # in it's own directory.
    megahit -1 {input.r1} -2 {input.r2} --min-contig-len 142 \
        --out-dir {wildcards.sample}_megahit \
        --out-prefix {wildcards.sample}
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.sample}_megahit/{wildcards.sample}.contigs.fa {output}
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.sample}_megahit
    '''

rule index_megahit:
    input: 'outputs/megahit/{sample}.contigs.fa'
    output: 'outputs/megahit/{sample}.contigs.fa.bwt'
    conda: 'bwa.yml'
    shell:'''
    bwa index {input}
    '''

rule map_to_megahit:
    input: 
        assembly='outputs/megahit/{sample}.contigs.fa',
        indx='outputs/megahit/{sample}.contigs.fa.bwt',
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz'
    output: 'outputs/map_to_megahit/{sample}.sam'
    conda:'bwa.yml'
    shell:'''
    bwa mem {input.assembly} {input.r1} {input.r2} > {output}
    '''

rule flagstat_megahit:
    input: 'outputs/map_to_megahit/{sample}.sam'
    output: 'outputs/map_to_megahit/{sample}.flagstat'
    conda:'samtools.yml'
    shell:'''
    samtools flagstat {input} > {output}
    '''
