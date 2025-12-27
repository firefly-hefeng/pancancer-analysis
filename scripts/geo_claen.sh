nohup python -u geo_cleaner.py >> geo_clean.log 2>&1 &
pid=$!
echo "$pid" | tee -a geo_clean.log
tail -f geo_clean.log