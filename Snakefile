# there's only 1 sample here: project_name.fastq.gz

# config:
#  ref:
#    sequence: /path/to/reference.fasta
#    genome_size: "20k"
#  hapdup:
#    sif: /varidata/research/projects/bbc/research/POSA_20240826_ONT_VBCS-1278/diploid/hapdup_0.12.sif
#    min_aligned_length: 5000
#  modules:
#    htslib: bbc2/htslib/htslib-1.17
#    minimap2: bbc2/minimap2/minimap2-2.26
#    samtools: bbc2/samtools/samtools-1.20

# import pandas as pd
configfile: "config/config.yaml"
project_name = config['project_name'] 
from pathlib import Path

rule all:
    input:
        # final stuff here
        # bam="results/align/{project_name}.bam"
        expand("results/subset/{project_name}_accepted.fastq.gz", project_name=project_name), #fill in here? SRCAP_PR001627_accepted
        expand("results/qc/all/NanoPlot-report-{project_name}-all.html", project_name=project_name),
        expand("results/qc/accepted/NanoPlot-report-{project_name}-accepted.html", project_name=project_name),
        expand("results/align/{project_name}.bam", project_name=project_name),
        expand("results/align/{project_name}_raw.bam", project_name=project_name),
        expand("results/align/{project_name}_subset.bam", project_name=project_name),
        expand("results/mosdepth/multiqc_report_{project_name}.mosdepth.html", project_name=project_name),
        expand("results/assemble/assembly-{project_name}.fasta", project_name=project_name),
        expand("results/phase/flye-{project_name}.bam", project_name=project_name),
        expand("results/phase/hapdup/hapdup_phased_1-{project_name}.fasta", project_name=project_name),
        expand("results/phase/hapdup/hapdup_phased_2-{project_name}.fasta", project_name=project_name),
        expand("results/deepvariant/{project_name}.vcf", project_name=project_name),
        expand("results/sniffles/{project_name}.vcf", project_name=project_name),
        expand("results/sniffles/{project_name}.PASS.called.nonrepeat.vcf", project_name=project_name),
        expand("results/sniffles/{project_name}.sniffles_igv_report.html", project_name=project_name)

rule subset:
    input:
        fq="raw_data/{project_name}.fastq.gz",
        csv="raw_data/adaptive_sampling_{project_name}.csv"
    output:
        # output files
        fq="results/subset/{project_name}_accepted.fastq.gz",
        csv="results/subset/{project_name}_accepted_reads.csv" #temp()
    envmodules:
        config['modules']['htslib'],
        config['modules']['seqtk'] 
    benchmark:
        "benchmarks/subset/{project_name}.txt"
    threads: 1
    resources: 
        mem_gb = 24,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        # find the lines in the csv that have "stop_receiving" + header
        ( head -n 1 {input.csv} && grep "sequence" $_ ) | sort | uniq > {output.csv}

        # subset the input fastq file
        # zgrep -A 3 -E $(echo -e $(tail -n +2 {output.csv} | awk -F, '{{print $1}}' | tr "\\n" "|" | sed "s/|$//")) {input.fq} | gzip -c > {output.fq}
        # zgrep fails bc argument list too long -- use seqtk instead

        # seqtk subseq Fastq.fq Header.lst > test.fq
        seqtk subseq {input.fq} <(tail -n +2 {output.csv} | awk -F, '{{print $1}}') | gzip -c > {output.fq}

        """


rule qc_all:
    input:
        fq="raw_data/{project_name}.fastq.gz"
    output:
        # output files
        all="results/qc/all/NanoPlot-report-{project_name}-all.html"
    benchmark:
        "benchmarks/qc/all-{project_name}.txt"
    conda:
        "nanoplot"
    threads: 8
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        NanoPlot \
            --threads {threads} \
            -o results/qc/all/ \
            --N50 \
            --fastq {input.fq}
        ln -sr results/qc/all/NanoPlot-report.html {output.all}
        """

rule qc_accepted:
    input:
        fq="results/subset/{project_name}_accepted.fastq.gz"
    output:
        accepted="results/qc/accepted/NanoPlot-report-{project_name}-accepted.html"
    benchmark:
        "benchmarks/qc/accepted-{project_name}.txt"
    conda:
        "nanoplot"
    threads: 8
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        NanoPlot \
            --threads {threads} \
            -o results/qc/accepted/ \
            --N50 \
            --fastq {input.fq}
        ln -sr results/qc/accepted/NanoPlot-report.html {output.accepted}
        """


rule align:
    input:
        "results/subset/{project_name}_accepted.fastq.gz"
    output:
        bam="results/align/{project_name}.bam"
    benchmark:
        "benchmarks/align/align-{project_name}.txt"
    params:
        ref_fasta=config["ref"]["sequence"]
    envmodules:
        config['modules']['minimap2'],
        config['modules']['samtools']
    threads: 128
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        minimap2 -t {threads} -x map-ont -a {params.ref_fasta} {input} | samtools sort -@ {threads} -o {output.bam} --write-index -
        """

rule align_raw:
    input:
        "raw_data/{project_name}.fastq.gz"
    output:
        bam="results/align/{project_name}_raw.bam"
    benchmark:
        "benchmarks/align/align_raw-{project_name}.txt"
    params:
        ref_fasta=config["ref"]["sequence"]
    envmodules:
        config['modules']['minimap2'],
        config['modules']['samtools']
    threads: 128
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        minimap2 -t {threads} -x map-ont -a {params.ref_fasta} {input} | samtools sort -@ {threads} -o {output.bam} --write-index -
        """

rule align_subset:
    input:
        "results/subset/{project_name}_accepted.fastq.gz"
    output:
        bam="results/align/{project_name}_subset.bam"
    benchmark:
        "benchmarks/align/align_subset-{project_name}.txt"
    params:
        ref_fasta_subset=config["ref"]["sequence_subset"]
    envmodules:
        config['modules']['minimap2'],
        config['modules']['samtools']
    threads: 128
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        minimap2 -t {threads} -x map-ont -a {params.ref_fasta_subset} {input} | samtools sort -@ {threads} -o {output.bam} --write-index -
        """

rule assemble:
    input:
        bam="results/align/{project_name}.bam",
        fq="results/subset/{project_name}_accepted.fastq.gz"
    output:
        assembly="results/assemble/assembly-{project_name}.fasta",
        assembly_info="results/assemble/assembly_info-{project_name}.txt"
    benchmark:
        "benchmarks/assemble/assemble-{project_name}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genome_size=config["ref"]["genome_size"]
    # envmodules:
    threads: 64
    resources: 
        mem_gb = 512,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        # until I make this BBC-accessible + an environmental module
        export PATH=$PATH:/varidata/research/projects/bbc/kaitlyn/tools/Flye/bin
        flye \
            --nano-raw {input.fq} \
            --genome-size {params.genome_size} \
            --out-dir $(dirname {output.assembly}) \
            --threads {threads}
        ln -sr $(dirname {output.assembly})/assembly.fasta {output.assembly}
        ln -sr $(dirname {output.assembly})/assembly_info.txt {output.assembly_info}
        """
    
rule phase:
    input:
        assembly="results/assemble/assembly-{project_name}.fasta",
        accepted_fastq="results/subset/{project_name}_accepted.fastq.gz",
    output:
        # phased="results/phase/phased-{project_name}.fasta",
        flye_bam="results/phase/flye-{project_name}.bam",
        phased_1_fasta="results/phase/hapdup/hapdup_phased_1-{project_name}.fasta",
        phased_2_fasta="results/phase/hapdup/hapdup_phased_2-{project_name}.fasta"
    benchmark:
        "benchmarks/phase/phase-{project_name}.txt"
    params:
        hapdup_sif=config["hapdup"]["sif"],
        hapdup_min_aligned_length=config["hapdup"]["min_aligned_length"]
    envmodules:
        config['modules']['minimap2'],
        config['modules']['samtools']
    threads: 120
    resources: 
        mem_gb = 960,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        """
        # cp {params.hapdup_sif} .

        #minimap, align to flye
        minimap2 -x map-ont -t {threads} -a {input.assembly} {input.accepted_fastq} | samtools sort -o {output.flye_bam} --write-index - 

        singularity exec --bind $PWD {params.hapdup_sif} hapdup \
            --assembly $PWD/{input.assembly} \
            --bam {output.flye_bam} \
            --out-dir $PWD/$(dirname {output.phased_2_fasta}) \
            -t {threads} \
            --rtype ont \
            --min-aligned-length {params.hapdup_min_aligned_length}

        tree -l -C 

        ln -sr $(dirname {output.phased_1_fasta})/hapdup_phased_1.fasta {output.phased_1_fasta}
        ln -sr $(dirname {output.phased_2_fasta})/hapdup_phased_2.fasta {output.phased_2_fasta}

        tree -l -C 

        sleep 30
        """

        # old dir, to bring back once testing is finished:
            # --out-dir $PWD/$(dirname {output.phased_2_fasta}) 
            # --out-dir $PWD/results/phase/hapdup/ \
        # ln -sr $(dirname {output.phased})/phased.fasta {output.phased}
        # ln -sr $(dirname {output.flye_bam})/flye.bam {output.flye_bam}
        # ln -sr {output.flye_bam} $(dirname {output.flye_bam})/flye.bam 

# rule depth:
#     input:
#         bam="results/align/{project_name}.bam",
#         ref_fasta=config["ref"]["sequence"]
#     output:
#         depth="results/depth/{project_name}.depth",
#         depth_multiqc="results/depth/multiqc_{project_name}.depth.html",
#     benchmark:
#         "benchmarks/depth/depth-{project_name}.txt"
#     envmodules:
#         config['modules']['samtools'],
#         config['modules']['multiqc']
#     threads: 16
#     resources: 
#         mem_gb = 64,
#         log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
#     shell:
#         # shell command
#         """
#         samtools depth -a -@ {threads} -o {output.depth} {input.bam}
#         samtools coverage 
#         multiqc {output.depth} \
#             --outdir results/depth/ \
#             --filename depth_{wildcards.project_name} \
#             --force
#         """

rule depth:
    input:
        bam="results/align/{project_name}.bam"
        # regions=config["ref"]["regions_depth_compute"]
    output:
        # mosdepth="results/mosdepth/{project_name}.mosdepth",
        fastq_summary="results/mosdepth/{project_name}.mosdepth.fast.mosdepth.summary.txt",
        # regions_folder="results/mosdepth/{project_name}.mosdepth.regions",
        multiqc_depth="results/mosdepth/multiqc_report_{project_name}.mosdepth.html"
    benchmark:
        "benchmarks/mosdepth/mosdepth-{project_name}.txt"
    params:
        fast_folder="results/mosdepth/{project_name}.mosdepth.fast"
    envmodules:
        config['modules']['mosdepth'],
        config['modules']['multiqc']
    threads: 16
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        # shell command
        # we've removed the following line: 
        # mosdepth --by construct.bed {output.regions_folder} {input.bam}
        """
        mosdepth -t {threads} -n --fast-mode {params.fast_folder} {input.bam}
        multiqc $(dirname {params.fast_folder}) \
            --filename {output.multiqc_depth} \
            --force
        """

rule sniffles:
    input:
        bam="results/align/{project_name}.bam",
        ref_fasta=config["ref"]["sequence"]
    output:
        vcf="results/sniffles/{project_name}.vcf"
    benchmark:
        "benchmarks/sniffles/sniffles-{project_name}.txt"
    params:
        repeats_bed=config["sniffles_repeats"]
    envmodules:
        config['modules']['sniffles']
    threads: 120
    resources: 
        mem_gb = 960,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        sniffles \
            -i {input.bam} \
            -v {output.vcf} \
            --reference {input.ref_fasta} \
            --tandem-repeats {params.repeats_bed} \
            --output-rnames \
            -t {threads}
        # python3 -m sniffles2_plot -i sniffles/SRCAP_PR001627.vcf -o sniffles
        """


def repeats_bed_path(ref_fasta):
    return f"results/sniffles/{Path(ref_fasta).stem}_repeats_merged.bed"
REPEATS_BED = repeats_bed_path(config["ref"]["sequence"])

rule prep_repeats_bed:
    input:
        ref_fasta=config["ref"]["sequence"]
    output:
        # repeats_bed="results/sniffles/{input.ref_fasta%.fa}_repeats_merged.bed"
        # repeats_bed=lambda wildcards, input:
        #     f"results/sniffles/{Path(input.ref_fasta).stem}_repeats_merged.bed"
        repeats_bed=REPEATS_BED #lambda wildcards, input: repeats_bed_path(input.ref_fasta)
    benchmark:
        "benchmarks/sniffles/prep_repeats_bed.txt"
    threads: 1
    resources: 
        mem_gb = 24,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    conda:
        "pbsv"
    shell:
        """
        # conda create -n pbsv bioconda::pbsv
        # conda activate pbsv
        /varidata/research/projects/bbc/kaitlyn/tools/pbsv/pbsv-2.11.0/annotations/findTandemRepeats \
            --merge \
            {input.ref_fasta} {output.repeats_bed}
        # per https://github.com/PacificBiosciences/pbsv/tree/master/annotations
        # and https://github.com/fritzsedlazeck/Sniffles/blob/master/annotations/source.txt

        """

rule filter_sniffles:
    input:
        vcf="results/sniffles/{project_name}.vcf",
        repeats_bed=lambda wildcards: repeats_bed_path(config["ref"]["sequence"])
    output:
        vcf_pass=temp("results/sniffles/{project_name}.PASS.vcf"),
        vcf_called=temp("results/sniffles/{project_name}.PASS.called.vcf"),
        positions_not_in_repeat_tsv=temp("results/sniffles/{project_name}.positions_not_in_repeat_bed.tsv"),
        vcf_no_repeats="results/sniffles/{project_name}.PASS.called.nonrepeat.vcf"
    benchmark:
        "benchmarks/sniffles/filter_sniffles-{project_name}.txt"
    envmodules:
        config['modules']['bedtools']
    threads: 1
    resources:
        mem_gb = 24,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        awk '$7 != "GT"' {input.vcf} > {output.vcf_pass}
        cat {output.vcf_pass} | grep -v "\./\." > {output.vcf_called}

        # get positions not in repeat.bed!
        bedtools intersect -v -a <(grep -v "^#" {output.vcf_called} | awk -vOFS="\t" '{{print $1, $2, $2}}' ) -b {input.repeats_bed} | cut -f 1-2 > {output.positions_not_in_repeat_tsv}

        (grep "^#" {output.vcf_called} && 
            while IFS= read -r line; do
                grep -w "$line" {output.vcf_called}; 
            done < {output.positions_not_in_repeat_tsv}) > {output.vcf_no_repeats}
        """

rule sniffles_igv_report_prep:
    input:
        vcf="results/sniffles/{project_name}.PASS.called.nonrepeat.vcf"
    output:
        bed="results/sniffles/{project_name}.PASS.called.nonrepeat.igv-reports.bed"
    benchmark:
        "benchmarks/sniffles/sniffles_igv_report_prep-{project_name}.txt"
    envmodules:
        config['modules']['bedtools']
    threads: 1
    resources:
        mem_gb = 24,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        # set up sniffles output bed for igv-reports
        awk -F'\t' 'BEGIN{{OFS="\t"}}
        $0 !~ /^#/ {{
        # Extract END from INFO
        if (match($8,/END=([^;]+)/,m)) {{
            end = m[1] + 1
        }} else {{
            end = $2 + 1
        }}

        # Extract END from INFO
        # end="."
        # if (match($8,/END=([^;]+)/,m)) end=m[1]
        
        # # Add +1 to END
        # end = end + 1

        # Build INFO: drop RNAMES, cut at STRAND, insert START=<POS> after SVLEN
        info=$8
        gsub(/RNAMES=[^;]+;?/,"",info)
        sub(/STRAND=.*$/,"",info)
        sub(/SVLEN=[^;]+;/,"&START=" $2 ";",info)
        gsub(/;+$/,"",info)

        # First three tab-separated fields
        printf "%s\t%s\t%s\t", $1,$2,end

        # Genotype/sample KEY:VALUE pairs (space-separated)
        nfmt = split($9, fmt, ":")
        nval = split($10, val, ":")
        for (i=1; i<=nfmt && i<=nval; i++) {{
            printf "%s:%s", fmt[i], val[i]
            if (i<nfmt && i<nval) printf " "
        }}

        # Space + INFO, then newline
        printf " %s\\n", info
        }}' {input.vcf} > {output.bed}
        """

rule sniffles_igv_report:
    input:
        bam="results/align/{project_name}.bam",
        bed="results/sniffles/{project_name}.PASS.called.nonrepeat.igv-reports.bed",
        fa=config["ref"]["sequence"],
        gtf=config["ref"]["gtf"]
    output:
        igv_report="results/sniffles/{project_name}.sniffles_igv_report.html"
    benchmark:
        "benchmarks/sniffles/sniffles_igv_report-{project_name}.txt"
    conda: "igv-reports"
    threads: 1
    resources:
        mem_gb = 24,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        create_report {input.bed} \
            --fasta {input.fa} \
            --title "Sniffles Structural Variants Report" \
            --flanking 100 \
            --window 23 \
            --roi {input.bed} \
            --standalone \
            --tracks {input.bam} {input.gtf} \
            --output {output.igv_report}
        """


rule deepvariant:
    input:
        bam="results/align/{project_name}.bam",
        ref_fasta=config["ref"]["sequence"]
    output:
        vcf="results/deepvariant/{project_name}.vcf",
        gvcf="results/deepvariant/{project_name}.gvcf",
        logs_dir=directory("results/deepvariant/logs-{project_name}")
    benchmark:
        "benchmarks/deepvariant/deepvariant-{project_name}.txt"
    envmodules:
        config['modules']['deepvariant']
    threads: 120
    resources:
        mem_gb = 960,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        # ln -s {input.ref_fasta} 
        # ln -s {input.ref_fasta}.fai 
        cp {input.ref_fasta} .
        cp {input.ref_fasta}.fai .

        echo $PWD $DEEPVAR
        tree -l -C $PWD
        singularity exec --bind $PWD/.. $DEEPVAR \
          /opt/deepvariant/bin/run_deepvariant \
          --model_type=ONT_R104 \
          --ref=$PWD/$(basename {input.ref_fasta}) \
          --reads={input.bam} \
          --output_vcf={output.vcf} \
          --output_gvcf={output.gvcf} \
          --num_shards={threads} \
          --vcf_stats_report=true \
          --logging_dir={output.logs_dir}
          
        # unlink $(basename {input.ref_fasta})
        # unlink $(basename {input.ref_fasta}.fai)

        # remove the copied files
        rm $(basename {input.ref_fasta})
        rm $(basename {input.ref_fasta}.fai)
        """

# singularity exec --bind /varidata/research/projects/bbc/research/LEMH_20250113_ONT_VBCS-1361/ $DEEPVAR \
#   /opt/deepvariant/bin/run_deepvariant \
#   --model_type=ONT_R104 \
#   --ref=../resources/fasta/mm39_chr7.fa \
#   --reads=results/align/SRCAP_PR001627.bam \
#   --output_vcf=deepvariant/output/ \
#   --output_gvcf=deepvariant/output/ \
#   --num_shards=$SLURM_CPUS_ON_NODE \
#   --vcf_stats_report=true \
#   --logging_dir=deepvariant/logs 

