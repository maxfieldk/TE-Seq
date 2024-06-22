Cytobands bed derived from:
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cytoBandIdeo.txt.gz
Telomere bed derived from:
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/gap.txt.gz
grep "telomere" gap.txt | cut -f2,3,4 | sort -k1,1V -k2,2n
cpgislands were derived from:
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cpgIslandExtUnmasked.txt.gz
cut -f 2,3,4,5,6,7,8,9,10 cpgIslandExtUnmasked.txt | sort -k1,1V -k2,2n -k3,3n > cpgislands.tsv