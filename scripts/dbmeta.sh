nohup python dbmeta_collect.py --all --output ./output > dbmeta.log 2>&1 &
PID=$!
echo "Rscript PID: $PID" >> dbmeta.log
tail -f dbmeta.log

nohup python cancer_infor.py > cancer_infor.log 2>&1 &


