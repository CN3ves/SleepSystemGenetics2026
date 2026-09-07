# SleepSystemGenetics2026
Workflow for the manuscript

Required packages:

```
module load miniforge3/25.3.0-3
conda activate snakemake
conda install python=3.13
conda install -c conda-forge -c bioconda deeptools
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
conda install bioconda::snakemake-executor-plugin-cluster-generic
conda list
```
```
required_packages <- c("BiocManager", "annotatr", "ATACseqQC", "biomaRt", "BSgenome.Mmusculus.UCSC.mm10", "ChIPseeker", "circlize", "clusterProfiler", "ComplexHeatmap", "csaw", "dendsort", "doParallel", "doSNOW", "dplyr", "EDASeq", "edgeR", "enrichplot", "foreach", "GenomicAlignments", "GenomicRanges", "ggplot2", "ggrepel", "GRaNIE", "limma", "optparse", "org.Mm.eg.db", "pheatmap", "preseqR", "qtl2", "RColorBrewer", "ReactomePA", "readr", "reshape2", "rjson", "Rsamtools", "Rsubread", "rtracklayer", "tidyr", "tidyverse", "TxDb.Mmusculus.UCSC.mm10.knownGene", "VennDiagram", "XML")

installed.packages <- installed.packages()[,"Package"]

for (package in required_packages) {
  if(! package %in% installed.packages){
    
    tryCatch({
      install.packages(package, lib='Rlibs')
    }, warning = function(cond) {
      BiocManager::install(package, lib='Rlibs')
    })
  }
}
```

Check DAG
```
snakemake -c 1 --dag  | dot -Tpng > workflow/dag.png; snakemake -c 1 --rulegraph  | dot -Tpng > workflow/rulegraph.png; snakemake -c 1 --filegraph  | dot -Tpng > workflow/filegraph.png; snakemake -c 1 -n; snakemake -n --debug-dag > snakemake_trace.txt 2>&1; snakemake --lint 2> lint; cat lint | grep 'Lints' | wc -l; cat lint | grep '*' | sort -u; rm lint
```

Run pipeline
```
snakemake --profile workflow/slurm_profile --use-conda --rerun-triggers code,input,mtime,params
```

```
echo "Installing Hmm-based IdeNtification of Transcription factor footprints" > {log}

#Need to downgrade everything
conda create -p ./rgt_env -c bioconda -c conda-forge python=3.11.7 # not enough space allocated
conda activate ./rgt_env

pip install RGT
pip install setuptools==69.5.1
pip install numpy==1.26.4
pip install scipy==1.12.0
  
cd ~/rgtdata
python setupGenomicData.py --mm10

conda deactivate 
```
```
rgt-hint  --version
HINT - Regulatory Analysis Toolbox (RGT) - v1.0.2
```

On remote servers, you might need to reinstall from time to time if you see a python version clash. run before
```
rm -r libs/python-3.11.7-hab00c5b_1_cpython* rgt_env/
```
If this happen, --rerun-triggers software-env will trigger a new analyses, hence why we used --rerun-triggers code,input,mtime,params (though ancient is used on inputs so it should render mtime useless)


<details>
  
  <summary>packages in environment at /users/USER/.conda/envs/snakemake:</summary>
  
Name | Version |  Build | Channel
|---|---|---|---|
_openmp_mutex |  4.5 | 20_gnu | conda-forge 
alabaster |  1.0.0 | pyhd8ed1ab_1 | conda-forge
amply |  0.1.7 | pyhd8ed1ab_0 | conda-forge
annotated-doc |  0.0.4 | pyhcf101f3_0 | conda-forge
annotated-types |  0.7.0 | pyhd8ed1ab_1 | conda-forge
argparse-dataclass | 2.0.0 | pyhd8ed1ab_1 | conda-forge
attrs |  26.1.0 |  pyhcf101f3_0 | conda-forge
babel |  2.18.0 |  pyhcf101f3_1 | conda-forge
backports.zstd | 1.6.0 |  py313h18e8e13_0 | conda-forge
boto3 |  1.43.52 | pypi_0 | pypi
botocore | 1.43.52 | pypi_0 | pypi
brotli | 1.2.0 | hed03a55_1 | conda-forge
brotli-bin | 1.2.0 | hb03c661_1 | conda-forge
brotli-python |  1.2.0 |  py313hf159716_1 | conda-forge
bzip2 |  1.0.8 | hda65f42_9 | conda-forge
c-ares | 1.34.6 |  hb03c661_0 | conda-forge
ca-certificates |  2026.6.17 | hbd8a1cb_0 | conda-forge
cairo |  1.18.4 |  he90730b_1 | conda-forge
certifi |  2026.6.17 | pyhd8ed1ab_0 | conda-forge
charset-normalizer | 3.4.7 | pyhd8ed1ab_0 | conda-forge
click |  8.4.2 | pypi_0 | pypi
coin-or-cbc |  2.10.13 | h4d16d09_1 | conda-forge
coin-or-cgl |  0.60.10 | hc46dffc_1 | conda-forge
coin-or-clp |  1.17.11 | hc03379b_1 | conda-forge
coin-or-osi |  0.108.12 |  hf4fecb4_1 | conda-forge
coin-or-utils |  2.11.13 | hc93afbd_1 | conda-forge
colorama | 0.4.6 | pyhd8ed1ab_1 | conda-forge
coloredlogs |  15.0.1 |  pyhd8ed1ab_4 | conda-forge
colormath2 | 3.0.3 | pypi_0 | pypi
conda-inject | 1.3.2 | pyhd8ed1ab_0 | conda-forge
configargparse | 1.7.5 | pyhcf101f3_0 | conda-forge
connection_pool |  0.0.3 | pyhd3deb0d_0 | conda-forge
contourpy |  1.3.3 |  py313hc8edb43_4 | conda-forge
cycler | 0.12.1 |  pyhcf101f3_2 | conda-forge
deeptools |  3.5.6 | pyhdfd78af_0 | bioconda
deeptoolsintervals | 0.1.9 |  py313hfeada96_12 | bioconda
docutils | 0.22.4 |  pyhd8ed1ab_0 | conda-forge
dpath |  2.2.0 | pyha770c72_1 | conda-forge
eido | 0.2.5 | pyhd8ed1ab_0 | conda-forge
font-ttf-dejavu-sans-mono | 2.37 |  hab24e00_0 | conda-forge
font-ttf-inconsolata | 3.000 | h77eed37_0 | conda-forge
font-ttf-source-code-pro | 2.038 | h77eed37_0 | conda-forge
font-ttf-ubuntu |  0.83 |  h77eed37_3 | conda-forge
fontconfig | 2.18.1 |  h27c8c51_0 | conda-forge
fonts-conda-ecosystem |  1 |  0 | conda-forge
fonts-conda-forge |  1 | hc364b38_1 | conda-forge
fonttools |  4.63.0 | py313h3dea7bd_0 | conda-forge
freetype | 2.14.3 |  ha770c72_0 | conda-forge
fribidi |  1.0.16 |  hb03c661_0 | conda-forge
gitdb |  4.0.12 |  pyhd8ed1ab_0 | conda-forge
gitpython |  3.1.50 |  pyhd8ed1ab_0 | conda-forge
graphite2 |  1.3.15 |  hecca717_0 | conda-forge
greenlet | 3.5.3 |  py313h5d5ffb9_0 | conda-forge
h2 | 4.3.0 | pyhcf101f3_0 | conda-forge
hpack |  4.2.0 | pyhd8ed1ab_0 | conda-forge
htseq |  2.1.2 | pypi_0 | pypi
humanfriendly |  10.0 |  pyh707e725_8 | conda-forge
humanize | 4.16.0 |  pypi_0 | pypi
hyperframe | 6.1.0 | pyhd8ed1ab_0 | conda-forge
icu |  78.3 |  h33c6efd_0 | conda-forge
idna | 3.18 |  pyhcf101f3_0 | conda-forge
imagesize |  2.0.0 | pyhd8ed1ab_0 | conda-forge
immutables | 0.21 | py313h07c4f96_2 | conda-forge
importlib-metadata | 9.0.0 | pyhcf101f3_0 | conda-forge
jinja2 | 3.1.6 | pyhcf101f3_1 | conda-forge
jmespath | 1.1.0 | pypi_0 | pypi
jsonschema | 4.26.0 |  pyhcf101f3_0 | conda-forge
jsonschema-specifications 2025.9.1 |  pyhcf101f3_0 | conda-forge
jupyter_core | 5.9.1 | pyhc90fa1f_0 | conda-forge
kaleido |  0.2.1 | pypi_0 | pypi
keyutils | 1.6.3 | hb9d3cd8_0 | conda-forge
kiwisolver | 1.5.0 |  py313hc8edb43_0 | conda-forge
krb5 | 1.22.2 |  hbde042b_1 | conda-forge
lcms2 |  2.19.1 |  h0c24ade_1 | conda-forge
ld_impl_linux-64 | 2.45.1 | default_hbd61a6d_102 | conda-forge
lerc | 4.1.0 | hdb68285_0 | conda-forge
libblas |  3.11.0 | 8_h4a7cf45_openblas | conda-forge
libbrotlicommon |  1.2.0 | hb03c661_1 | conda-forge
libbrotlidec | 1.2.0 | hb03c661_1 | conda-forge
libbrotlienc | 1.2.0 | hb03c661_1 | conda-forge
libcblas | 3.11.0 | 8_h0358290_openblas | conda-forge
libcurl |  8.21.0 |  hae6b9f4_2 | conda-forge
libdeflate | 1.25 |  h17f619e_0 | conda-forge
libedit |  3.1.20250104 | pl5321h7949ede_0 | conda-forge
libev |  4.33 |  hd590300_2 | conda-forge
libexpat | 2.8.1 | hecca717_1 | conda-forge
libffi | 3.5.2 | h3435931_0 | conda-forge
libfreetype |  2.14.3 |  ha770c72_0 | conda-forge
libfreetype6 | 2.14.3 |  h73754d4_0 | conda-forge
libgcc | 15.2.0 | he0feb66_19 | conda-forge
libgcc-ng |  15.2.0 | h69a702a_19 | conda-forge
libgfortran |  15.2.0 | h69a702a_19 | conda-forge
libgfortran5 | 15.2.0 | h68bc16d_19 | conda-forge
libglib |  2.88.2 |  h0d30a3d_0 | conda-forge
libgomp |  15.2.0 | he0feb66_19 | conda-forge
libharfbuzz |  14.2.1 |  h17a8019_1 | conda-forge
libiconv | 1.18 |  h3b78370_2 | conda-forge
libjpeg-turbo |  3.1.4.1 | hb03c661_0 | conda-forge
liblapack |  3.11.0 | 8_h47877c9_openblas | conda-forge
liblapacke | 3.11.0 | 8_h6ae95b6_openblas | conda-forge
liblzma |  5.8.3 | hb03c661_0 | conda-forge
libmpdec | 4.0.0 | hb03c661_1 | conda-forge
libnghttp2 | 1.68.1 |  h877daf1_0 | conda-forge
libopenblas |  0.3.33 | pthreads_h94d23a6_0 | conda-forge
libpng | 1.6.58 |  h421ea60_0 | conda-forge
libpsl | 0.22.0 |  hd9031aa_0 | conda-forge
libraqm |  0.10.5 |  h6406941_1 | conda-forge
libsqlite |  3.53.3 |  h0c1763c_0 | conda-forge
libssh2 |  1.11.1 |  hcf80075_0 | conda-forge
libstdcxx |  15.2.0 | h934c35e_19 | conda-forge
libstdcxx-ng | 15.2.0 | hdf11a46_19 | conda-forge
libtiff |  4.7.2 | h9d88235_0 | conda-forge
libuuid |  2.42.2 |  h5347b49_0 | conda-forge
libwebp-base | 1.6.0 | hd42ef1d_0 | conda-forge
libxcb | 1.17.0 |  h8a09558_0 | conda-forge
libzlib |  1.3.2 | h25fd6f3_2 | conda-forge
logmuse |  0.3.0 | pyhcf101f3_0 | conda-forge
markdown | 3.10.2 |  pypi_0 | pypi
markdown-it-py | 4.2.0 | pyhd8ed1ab_0 | conda-forge
markupsafe | 3.0.3 |  py313h3dea7bd_1 | conda-forge
matplotlib-base |  3.11.0 | py313h6c470cf_1 | conda-forge
mdurl |  0.1.2 | pyhd8ed1ab_1 | conda-forge
multiqc |  1.35 |  pypi_0 | pypi
munkres |  1.1.4 | pyhd8ed1ab_1 | conda-forge
narwhals | 2.23.0 |  pyhcf101f3_0 | conda-forge
natsort |  8.4.0 | pypi_0 | pypi
nbformat | 5.10.4 |  pyhd8ed1ab_1 | conda-forge
ncurses |  6.6 | hdb14827_0 | conda-forge
networkx | 3.6.1 | pypi_0 | pypi
numpy |  2.5.1 |  py313hf6604e3_0 | conda-forge
numpydoc | 1.10.0 |  pyhcf101f3_0 | conda-forge
openjpeg | 2.5.4 | h55fea9a_0 | conda-forge
openssl |  3.6.3 | h35e630c_0 | conda-forge
packaging |  25.0 |  pyh29332c3_1 | conda-forge
pandas | 2.3.3 |  py313h08cd8bf_2 | conda-forge
pcre2 |  10.47 | haa7fec5_0 | conda-forge
pephubclient | 0.4.4 | pyhd8ed1ab_1 | conda-forge
peppy |  0.40.8 |  pyhd8ed1ab_0 | conda-forge
pillow | 12.3.0 | py313h80991f8_0 | conda-forge
pip |  26.1.2 |  pyh145f28c_0 | conda-forge
pixman | 0.46.4 |  h54a6638_1 | conda-forge
platformdirs | 4.10.0 |  pyhcf101f3_0 | conda-forge
plotly | 6.8.0 | pyhd8ed1ab_0 | conda-forge
polars | 1.43.0 |  pypi_0 | pypi
polars-runtime-32 |  1.43.0 |  pypi_0 | pypi
polars-runtime-compat |  1.43.0 |  pypi_0 | pypi
psutil | 7.2.2 |  py313h54dd161_0 | conda-forge
pthread-stubs |  0.4 |  hb9d3cd8_1002 | conda-forge
pulp | 2.8.0 |  py313hf1034c9_3 | conda-forge
py2bit | 1.0.1 |  py313hd978853_0 | bioconda
pyarrow |  25.0.0 |  pypi_0 | pypi
pybigwig | 0.3.25 | py313h759994b_1 | bioconda
pydantic | 2.13.4 |  pyhcf101f3_0 | conda-forge
pydantic-core |  2.46.4 | py313h843e2db_0 | conda-forge
pygments | 2.20.0 |  pyhd8ed1ab_0 | conda-forge
pyparsing |  3.3.2 | pyhcf101f3_0 | conda-forge
pysam |  0.24.0 | py313h4b224ce_1 | bioconda
pysocks |  1.7.1 | pyha55dd90_7 | conda-forge
python | 3.13.14 |  h6add32d_100_cp313 | conda-forge
python-dateutil |  2.9.0.post0 | pyhe01879c_2 | conda-forge
python-dotenv |  1.2.2 | pypi_0 | pypi
python-fastjsonschema |  2.21.2 |  pyhe01879c_0 | conda-forge
python-tzdata |  2026.2 |  pyhd8ed1ab_0 | conda-forge
python_abi | 3.13 | 8_cp313 | conda-forge
pytz | 2026.2 |  pyhcf101f3_0 | conda-forge
pyyaml | 6.0.3 |  py313h3dea7bd_1 | conda-forge
qhull |  2020.2 |  h434a139_5 | conda-forge
readline | 8.3 | h853b02a_0 | conda-forge
referencing |  0.37.0 |  pyhcf101f3_0 | conda-forge
regex |  2026.7.19 | pypi_0 | pypi
requests | 2.34.2 |  pyhcf101f3_0 | conda-forge
rich | 15.0.0 |  pyhcf101f3_0 | conda-forge
rich-click | 1.9.8 | pypi_0 | pypi
roman-numerals | 4.1.0 | pyhd8ed1ab_0 | conda-forge
rpds-py |  2026.6.3 | py313hafbe609_0 | conda-forge
s3transfer | 0.19.1 |  pypi_0 | pypi
scipy |  1.18.0 | py313h4b8bb8b_0 | conda-forge
shellingham |  1.5.4 | pyhd8ed1ab_2 | conda-forge
six |  1.17.0 |  pyhe01879c_1 | conda-forge
slack-sdk |  3.43.0 |  pyhcf101f3_0 | conda-forge
slack_sdk |  3.43.0 |  pyh9dfb50f_0 | conda-forge
smart_open | 7.7.1 | pyhcf101f3_0 | conda-forge
smmap |  5.0.3 | pyhcf101f3_1 | conda-forge
snakemake |  9.23.1 |  hdfd78af_1 | bioconda
snakemake-executor-plugin-cluster-generic | 1.0.9 | pyhdfd78af_0 | bioconda
snakemake-interface-common  |  1.23.0 |  pyhdfd78af_1 | bioconda
snakemake-interface-executor-plugins  | 9.4.0 | pyh84498cf_0 | bioconda
snakemake-interface-logger-plugins | 2.1.0 | pyhdfd78af_0 | bioconda
snakemake-interface-report-plugins | 1.3.0 | pyhd4c3c12_0 | bioconda
snakemake-interface-scheduler-plugins | 2.0.2 | pyhd4c3c12_0 | bioconda
snakemake-interface-storage-plugins | 4.4.1 | pyh84498cf_0 | bioconda
snakemake-minimal |  9.23.1 |  pyhdfd78af_1 | bioconda
snowballstemmer |  3.1.1 | pyhd8ed1ab_0 | conda-forge
spectra |  0.1.0 | pypi_0 | pypi
sphinx | 9.1.0 | pyhd8ed1ab_0 | conda-forge
sphinxcontrib-applehelp |  2.0.0 | pyhd8ed1ab_1 | conda-forge
sphinxcontrib-devhelp |  2.0.0 | pyhd8ed1ab_1 | conda-forge
sphinxcontrib-htmlhelp | 2.1.0 | pyhd8ed1ab_1 | conda-forge
sphinxcontrib-jsmath | 1.0.1 | pyhd8ed1ab_1 | conda-forge
sphinxcontrib-qthelp | 2.0.0 | pyhd8ed1ab_1 | conda-forge
sphinxcontrib-serializinghtml | 2.0.0 | pyhd8ed1ab_0 | conda-forge
sqlalchemy | 2.0.51 | py313h54dd161_0 | conda-forge
sqlmodel | 0.0.37 |  pyhcf101f3_0 | conda-forge
tabulate | 0.10.0 |  pyhcf101f3_0 | conda-forge
tenacity | 9.1.4 | pyhcf101f3_0 | conda-forge
throttler |  1.2.2 | pyhd8ed1ab_0 | conda-forge
tiktoken | 0.13.0 |  pypi_0 | pypi
tk | 8.6.13 | noxft_h366c992_103 | conda-forge
toml | 0.10.2 |  pyhcf101f3_3 | conda-forge
tomli |  2.4.1 | pyhcf101f3_0 | conda-forge
tqdm | 4.69.0 |  pypi_0 | pypi
traitlets |  5.15.1 |  pyhcf101f3_0 | conda-forge
typeguard |  4.5.2 | pypi_0 | pypi
typer |  0.26.8 |  pyhcf101f3_0 | conda-forge
typing-extensions |  4.16.0 |  h69aa097_0 | conda-forge
typing-inspection |  0.4.2 | pyhcf101f3_2 | conda-forge
typing_extensions |  4.16.0 |  pyhcf101f3_0 | conda-forge
tzdata | 2025c | hc9c84f9_1 | conda-forge
ubiquerg | 0.9.3 | pyhd8ed1ab_0 | conda-forge
urllib3 |  2.7.0 | pyhd8ed1ab_0 | conda-forge
wrapt |  1.17.3 | py313h07c4f96_1 | conda-forge
xorg-libice |  1.1.2 | hb9d3cd8_0 | conda-forge
xorg-libsm | 1.2.6 | he73a12e_0 | conda-forge
xorg-libx11 |  1.8.13 |  he1eb515_0 | conda-forge
xorg-libxau |  1.0.12 |  hb03c661_1 | conda-forge
xorg-libxdmcp |  1.1.5 | hb03c661_1 | conda-forge
xorg-libxext | 1.3.7 | hb03c661_0 | conda-forge
xorg-libxrender |  0.9.12 |  hb9d3cd8_0 | conda-forge
yaml | 0.2.5 | h280c20c_3 | conda-forge
yte |  1.9.4 | pyhd8ed1ab_0 | conda-forge
zipp | 4.1.0 | pyhcf101f3_0 | conda-forge
zlib | 1.3.2 | h25fd6f3_2 | conda-forge
zlib-ng |  2.3.3 | hceb46e0_1 | conda-forge
zstd | 1.5.7 | hb78ec9c_6 | conda-forge
</details>

<details>
  
  <summary>R packages loaded:</summary>

Package | Version
|---|---|
foreach | 1.5.2 
doParallel | 1.0.17 
doSNOW | 1.0.20 
rtracklayer | 1.70.1 
colorspace | 2.1-3
csaw | 1.44.0 
tidyverse | 2.0.0 
GenomicRanges | 1.62.1 
optparse | 1.8.2 
ChIPseeker | 1.46.1 
ggplot2 | 4.0.3 
dplyr | 1.2.1 
BSgenome.Mmusculus.UCSC.mm10 | 1.4.3 
TxDb.Mmusculus.UCSC.mm10.knownGene | 3.10.0 
ggrepel | 0.9.8 
annotatr | 1.36.0 
ComplexHeatmap | 2.26.1 
circlize | 0.4.18 
org.Mm.eg.db | 3.22.0 
clusterProfiler | 4.18.4 
enrichplot | 1.30.5 
ReactomePA | 1.54.0 
RColorBrewer | 1.1-3 
dendsort | 0.3.4 
edgeR | 4.8.2 
VennDiagram | 1.8.2 
readr | 2.2.0 
GRaNIE | 1.14.0 
limma | 3.66.0 
biomaRt | 2.66.2 
rjson | 0.2.23 
EDASeq | 2.44.0 
ATACseqQC | 1.34.0 
Rsamtools | 2.26.0 
openxlsx | 4.2.8.1
preseqR | 4.0.0 
GenomicAlignments | 1.46.0 
tidyr | 1.3.2 
reshape2 | 1.4.5 
XML | 3.99-0.24 
Rsubread | 2.24.0 
pheatmap | 1.0.13 
qgraph | 2.3.3
igraph | 2.3.3
qtl2 | 0.46 
svglite | 2.2.2
</details>

```
x <- read.csv('results/12-BXD_footprints/plots/footprint_analysis.csv')
x <-unique(unlist(strsplit(x[x$padj<0.05,1],':')))
x <-  gsub('.*\\.','',gsub('\\(.*','',x))
cat(paste(sort(unique(toupper(x))), collapse='\n'))
```
https://maayanlab.cloud/kea3/

place in results/12-BXD_footprints/KEA/Integrated scaled rank.tsv
place in results/12-BXD_footprints/KEA/Mean rank.tsv

some rules gather results from other rules so they might need to be run sequentially. the pipeline was testes runnign one rule at a time, not the entire workflow