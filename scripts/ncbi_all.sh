nohup python -u ncbi_all_collector.py >> ncbi_all.log 2>&1 &
pid=$!
echo "$pid" | tee -a ncbi_all.log
tail -f ncbi_all.log