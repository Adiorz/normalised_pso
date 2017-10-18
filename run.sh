#!/bin/bash

num_modes=5
num_iter=50
num_particles=50

c=$1

data_file="data/04_mod_A20_FRF_displacement.csv"

pdfs_dir="pdfs_${num_modes}_modes_04_mod_A20^2"
logs_dir="logs_${num_modes}_modes_04_mod_A20^2"
pdf_file="channel_${c}.pdf"
log_file="channel_${c}.log"

python_path="PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages"

echo "$pdfs_dir"
echo "$logs_dir"
mkdir -p "$pdfs_dir"
mkdir -p "$logs_dir"

#rm -f "$logs_dir/04_mod_A20_FRF_log_channel_$c.log"
#rm -f "$logs_dir/04_mod_A20_FRF_log_channel_$c.log_main"

#debug
#gdb -ex run --args ./main "$data_file" "$logs_dir//$log_file" "$c" 1001 "$num_modes" "$num_iter" "$num_particles"

#release
#./main "$data_file" "$logs_dir/$log_file" "$c" 1001 "$num_modes" "$num_iter" "$num_particles"

#main workers results
#echo "$python_path ./plot.py $logs_dir/${log_file}_main $pdfs_dir/$pdf_file"
#PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs_dir/${log_file}_main" "$pdfs_dir/$pdf_file"

#helpers results
#echo "$python_path ./plot.py $logs_dir/$log_file $pdfs_dir/pdf_file"
#PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs_dir/$log_file" "$pdfs_dir/pdf_file"

#for c in {@..16}
for c in {0..0}
do
    echo "$c"
    pdf_file="channel_${c}.pdf"
    log_file="channel_${c}.log"
    echo "$log_file, $pdf_file"
    ./main "$data_file" "$logs_dir/$log_file" "$c" 1001 "$num_modes" "$num_iter" "$num_particles"
    PYTHONPATH=/home/adiorz/anaconda3/lib/python3.5/site-packages ./plot.py "$logs_dir/$log_file" "$pdfs_dir/$pdf_file"
done
