#!/bin/bash
nohup pypy /work-z/user/guoh/tools/findVar.py -f Tumor_loci_Colon-Lung_cancer.txt -v SSCS_pyvacheck/ -o sscs_pyvarcheck.txt >findvar.out&
nohup pypy /work-z/user/guoh/tools/findVar.py -f Tumor_loci_Colon-Lung_cancer.txt -v pyvacheck/ -o pyvarcheck.txt >findvar.out&
