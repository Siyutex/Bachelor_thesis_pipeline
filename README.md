# Bachelor\_thesis\_pipeline

This repository contains the code I write for my bachelor thesis project involving boolean network inference on scRNAseq data and subsequent attractor analysis

IMPORTANT INFO: 
- the  "helper\_functions" library is a personal library not available on pypi. If you want to use the scripts as is, you can download it with:
	pip install --index-url https://test.pypi.org/simple/ JL\_helper\_functions

-there is a bug in arboreto 0.1.6 that yields the error "must supply at least on delayed object"; to fix go into arboreto.core.create\_graph and comment out line the line that says "all\_meta\_df = from\_delayed(delayed\_meta\_dfs, meta=\_META\_SCHEMA)"

