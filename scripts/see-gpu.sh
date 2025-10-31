#!/usr/bin/env bash
echo "========== CPU =========="
lscpu | grep -E 'Model name|Socket|Core|Thread'
echo
echo "========== Memory =========="
free -h
echo
echo "========== GPU =========="
nvidia-smi --query-gpu=name,driver_version,memory.total,memory.used,temperature.gpu,utilization.gpu \
           --format=csv,noheader,nounits
echo
echo "========== Disk =========="
df -h / | tail -1
echo
echo "========== CUDA / PyTorch 可见性 =========="
python -c "import torch, sys; print('Python', sys.version.split()[0]); \
           print('PyTorch', torch.__version__); \
           print('CUDA available:', torch.cuda.is_available(), \
                 '| device count:', torch.cuda.device_count(), \
                 '| current:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'None')"