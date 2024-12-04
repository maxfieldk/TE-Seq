coords=$(awk '{print $1":"$2"-"$3}' <<< $(grep "L1HS" /users/mkelsey/data/Nanopore/alz/RTE/ldna/Rintermediates/flRTEpromoter.tsv))
rte_name=$(awk '{print $6}' <<< $(grep "L1HS" /users/mkelsey/data/Nanopore/alz/RTE/ldna/Rintermediates/flRTEpromoter.tsv))
mapfile -t coords_array <<< "$coords"
mapfile -t rte_name_array <<< "$rte_name"
for (( i=0; i<${#coords_array[@]}; i++ ))
do
coord=${coords_array[$i]}
rte_name=${rte_name_array[$i]}
coordfnsafe=$(echo $coord | sed 's/:/-/')
/users/mkelsey/data/tools/wgbs_tools/wgbstools pat_fig --nanopore --shuffle --line_width 0.5 --col_wrap 1 --font_size 3 --meth_color "beige" --unmeth_color "black"  -o $(dirname ldna/results/wgbs/plots/l1hs/by_sample/l1hs_fl_by_sample.outfile)/${rte_name}.pdf -r $coord --title $rte_name --strict --genome hs1 $(dirname $(dirname ldna/results/wgbs/analysis_default/merged/merged.outfile))/*.pat.gz
done
touch {output.outfile}