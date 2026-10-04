# Automatic Loupe alignment
```shell
cd /path/to/Spa-PR/PR/P3/0_PR_03-loupe-Spaceranger-results
export PATH=/path/to/biosofts/spaceranger-2.0.0:$PATH
spaceranger count --id="PR_P3_Spaceranger" \
                 --description="PR (FFPE)" \
                 --transcriptome=/path/to/ref/GRCh38-ref-for-space-ranger-latest/refdata-gex-GRCh38-2020-A \
                 --probe-set=/path/to/ref/GRCh38-ref-for-space-ranger-latest/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \
                 --fastqs=/path/to/Spa-PR/PR/P3/data/ \
                 --sample=PRAD_385-1 \
                 --image=/path/to/Spa-PR/PR/P3/data/PRAD_385.jpg \
                 --slide=V11N29-070 \
                 --area=D1 \
                 --reorient-images=true \
                 --localcores=16 \
                 --localmem=64
```
# Manual Loupe alignment
```shell
cd /path/to/Spa-PR/PR/P3/0_PR_03-loupe-Spaceranger-results
export PATH=/path/to/biosofts/spaceranger-2.0.0:$PATH
spaceranger count --id="PR_P3_Spaceranger" \
                 --description="PR (FFPE)" \
                 --transcriptome=/path/to/ref/GRCh38-ref-for-space-ranger-latest/refdata-gex-GRCh38-2020-A \
                 --probe-set=/path/to/ref/GRCh38-ref-for-space-ranger-latest/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \
                 --fastqs=/path/to/Spa-PR/PR/P3/data/ \
                 --sample=PRAD_385-1 \
                 --image=/path/to/Spa-PR/PR/P3/data/PRAD_385.jpg \
                 --slide=V11N29-070 \
                 --area=D1 \
                 --localcores=16 \
                 --localmem=64 \
                 --loupe-alignment=/path/to/Spa-PR/PR/P3/data/V11N29-070-D1.json
```
