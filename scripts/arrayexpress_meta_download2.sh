nohup python -u arrayexpress_meta_download2.py >> arrayexpress_meta_download2.log 2>&1 &
pid=$!
echo "$pid" | tee -a arrayexpress_meta_download2.log
tail -f arrayexpress_meta_download2.log