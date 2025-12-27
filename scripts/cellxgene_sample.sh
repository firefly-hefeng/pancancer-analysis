nohup python -u cellxgene_sample.py >> cellxgene_sample.log 2>&1 &
pid=$!
echo "$pid" | tee -a cellxgene_sample.log
tail -f cellxgene_sample.log