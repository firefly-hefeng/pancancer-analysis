nohup python -u geo_reviewer.py >> geo_review.log 2>&1 &
pid=$!
echo "$pid" | tee -a geo_review.log
tail -f geo_review.log