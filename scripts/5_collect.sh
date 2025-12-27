nohup python -u dryad_collector.py >> dryad.log 2>&1 &
pid=$!
echo "$pid" | tee -a dryad.log
nohup python -u zenodo_collector.py >> zenodo.log 2>&1 &
pid=$!
echo "$pid" | tee -a zenodo.log
nohup python -u figshare_collector.py >> figshare.log 2>&1 &
pid=$!
echo "$pid" | tee -a figshare.log
nohup python -u broad_scp_collector.py >> broad_scp.log 2>&1 &
pid=$!
echo "$pid" | tee -a braod_scp.log
nohup python -u cellxgene_collector.py >> cellxgene.log 2>&1 &
pid=$!
echo "$pid" | tee -a cellxgene.log

