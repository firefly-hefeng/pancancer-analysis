nohup python -u meta1_cleaner.py >> meta1_clean.log 2>&1 &
pid=$!
echo "$pid" | tee -a meta1_clean.log
tail -f meta1_clean.log