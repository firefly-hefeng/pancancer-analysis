nohup python -u arrayexpress_collector.py >> arrayexpress.log 2>&1 &
pid=$!
echo "$pid" | tee -a arrayexpress.log
tail -f arrayexpress.log