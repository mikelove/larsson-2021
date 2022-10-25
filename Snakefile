CELLS = range(1,384)

SALMON = "/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon"

ANNO = "/proj/milovelab/love/proj/ase/osteoblast-quant/diploid_txomes/indices/CAST_EiJ"

rule all:
    input: expand("quants/cell{cell}/quant.sf", cell=CELLS)

rule salmon_quant:
    input:
        "fastq/{cell}"
    output:
        "quants/{cell}/quant.sf"
    params:
        dir = "quants/{cell}"
    shell:
        "{SALMON} quant -i {ANNO} -l U -p 12 --numBootstraps 30 -o {params.dir} -r {input}/*.fastq.gz"
