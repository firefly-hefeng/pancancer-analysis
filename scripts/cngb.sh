nohup python -u cngb_collector.py >> cngb.log 2>&1 &
pid=$!
echo "$pid" | tee -a cngb.log
tail -f cngb.log