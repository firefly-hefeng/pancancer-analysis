nohup python -u meta0_cleaner.py >> meta0_clean.log 2>&1 &
pid=$!
echo "$pid" | tee -a meta0_clean.log
tail -f meta0_clean.log