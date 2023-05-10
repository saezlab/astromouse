
from snakemake.remote import AUTO

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()


rule deconvolution:
    input:
        HTTP.remote('data.mendeley.com/api/datasets/fjxrcbh672/draft/files/7c5f83e8-8fdd-4818-822a-ad514659676f?a=69394d54-235c-436e-be60-520cd2899517', keep_local=False)
    output:
        'data/original/ST/ST_brain_deconvoluted.rds'
    shell:
        "(test -d data/ || mkdir data) && "
        "(test -d data/original/ || mkdir data/original) && "
        "(test -d data/original/ST/ || mkdir data/original/ST) && "
        "(test -d temp/ || mkdir temp) && "
        "unzip {input} deconvolution_results/sobject_brain_deconv_220609.Rds -d temp && "
        "mv temp/deconvolution_results/sobject_brain_deconv_220609.Rds {output} && "
        "rm -r temp/"

checkpoint visium_data:
    input:
        HTTP.remote("data.mendeley.com/api/datasets/fjxrcbh672/draft/files/ab629966-ac1a-4b9e-8371-a739a14a859b?a=69394d54-235c-436e-be60-520cd2899517", keep_local=False)
    output:
        directory('data/original/ST/visium_data_brain')
    shell:
        "(test -d data/ || mkdir data) && "
        "(test -d data/original/ || mkdir data/original) && "
        "(test -d data/original/ST/ || mkdir data/original/ST) && "
        "(test -d temp/ || mkdir temp) && "
        "unzip {input} 'visium_data/*' -d temp && "
        "mv temp/visium_data {output} && "
        "rm -r temp/"

rule multiome_data:
    output:
        'data/original/MO/MO_brain_annotated.RData'
    shell:
        "(test -d data/ || mkdir data) && "
        "(test -d data/original/ || mkdir data/original) && "
        "(test -d data/original/MO/ || mkdir data/original/MO) && "
        "(test -d temp/ || mkdir temp) && "
        "curl -L 'https://data.mendeley.com/api/datasets/fjxrcbh672/draft/files/43ae9888-e8a4-4d96-8098-48d7506b549b?a=69394d54-235c-436e-be60-520cd2899517' -o temp/MO.zip && "
        "unzip temp/MO.zip multiomics_final_objects_metadata/seurat_data_brain_cc_2022-12-22.RData -d temp && "
        "mv temp/multiomics_final_objects_metadata/seurat_data_brain_cc_2022-12-22.RData {output} "
        "rm -r temp/"

rule ST_annotated:
    input:
        HTTP.remote('https://data.mendeley.com/api/datasets/fjxrcbh672/draft/files/122cc911-9eb7-481e-83da-1b3f70cb88d3?a=69394d54-235c-436e-be60-520cd2899517', keep_local=False)
    output:
        'data/original/ST/ST_brain_annotated.rds'
    shell:
        "(test -d data/ || mkdir data) && "
        "(test -d data/original/ || mkdir data/original) && "
        "(test -d data/original/ST/ || mkdir data/original/ST) && "
        "mv {input} {output}"