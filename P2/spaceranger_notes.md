# Automatic Loupe alignment
```shell
cd /path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results
export PATH=/path/to/biosofts/spaceranger-2.0.0:$PATH
spaceranger count --id="PR_P2_Spaceranger" \
                 --description="PR (FFPE)" \
                 --transcriptome=/path/to/ref/GRCh38-ref-for-space-ranger-latest/refdata-gex-GRCh38-2020-A \
                 --probe-set=/path/to/ref/GRCh38-ref-for-space-ranger-latest/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \
                 --fastqs=/path/to/Spa-PR/PR/P2/data/ \
                 --sample=PRAD_3007-1,PRAD_3007-2 \
                 --image=/path/to/Spa-PR/PR/P2/data/PRAD_3007.jpg \
                 --slide=V11N29-073 \
                 --area=C1 \
                 --reorient-images=true \
                 --localcores=16 \
                 --localmem=64
```
# Manual Loupe alignment
```shell
cd /path/to/Spa-PR/PR/P2/0_PR_02-loupe-Spaceranger-results
export PATH=/path/to/biosofts/spaceranger-2.0.0:$PATH
spaceranger count --id="PR_P2_Spaceranger" \
                 --description="PR (FFPE)" \
                 --transcriptome=/path/to/ref/GRCh38-ref-for-space-ranger-latest/refdata-gex-GRCh38-2020-A \
                 --probe-set=/path/to/ref/GRCh38-ref-for-space-ranger-latest/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \
                 --fastqs=/path/to/Spa-PR/PR/P2/data/ \
                 --sample=PRAD_3007-1,PRAD_3007-2 \
                 --image=/path/to/Spa-PR/PR/P2/data/PRAD_3007.jpg \
                 --slide=V11N29-073 \
                 --area=C1 \
                 --localcores=16 \
                 --localmem=64 \
                 --loupe-alignment=/path/to/Spa-PR/PR/P2/data/V11N29-073-C1.json
```
