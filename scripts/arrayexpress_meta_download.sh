nohup python -u arrayexpress_meta_download.py >> arrayexpress_meta_download.log 2>&1 &
pid=$!
echo "$pid" | tee -a arrayexpress_meta_download.log
tail -f arrayexpress_meta_download.log