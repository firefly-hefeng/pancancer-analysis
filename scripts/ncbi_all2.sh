nohup python -u ncbi_all_2.py >> ncbi_all2.log 2>&1 &
pid=$!
echo "$pid" | tee -a ncbi_all2.log
tail -f ncbi_all2.log