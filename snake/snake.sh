#!/usr/bin/bash
snakemake -s CRISPRscreen_snake.py --cluster "qsub -l nodes=1:ppn={threads} -q all -V -S /bin/bash -p 1000" --jobs 432 --keep-going --rerun-incomplete -R --max-jobs-per-second 432 -p

