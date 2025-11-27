import subprocess
import shutil

def get_gpu_info():
    """
    Detects available GPUs using nvidia-smi.
    Returns a list of dictionaries with 'name', 'total_memory_mb', 'free_memory_mb'.
    """
    if not shutil.which("nvidia-smi"):
        return []

    try:
        # Run nvidia-smi to get GPU details
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total,memory.free", "--format=csv,noheader,nounits"],
            capture_output=True,
            text=True,
            check=True
        )
        
        gpus = []
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if not line:
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 3:
                gpus.append({
                    "name": parts[0],
                    "total_memory_mb": int(parts[1]),
                    "free_memory_mb": int(parts[2])
                })
        return gpus
    except Exception as e:
        print(f"Error detecting GPUs: {e}")
        return []

def estimate_max_workers(gpus, model_size_gb=10):
    """
    Estimates the maximum number of concurrent workers based on free VRAM.
    Assumes model_size_gb is the VRAM required for one instance/stream of the model.
    """
    if not gpus:
        return 1 # Fallback for CPU
        
    total_free_vram_gb = sum(gpu['free_memory_mb'] for gpu in gpus) / 1024
    
    # Simple heuristic: Total Free VRAM / Model Size
    # We add a small buffer or floor it.
    estimated_workers = int(total_free_vram_gb // model_size_gb)
    
    return max(1, estimated_workers)
