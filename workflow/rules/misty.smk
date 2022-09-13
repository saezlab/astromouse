rule get_coords:
    input:
        data = 'data/working/ST/{tissue}_wImages.h5ad'
    output:
        coords = 'data/working/ST/Misty/{tissue}_coordinates.csv'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/misty/extract_coordinates.py"

