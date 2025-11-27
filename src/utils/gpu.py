import subprocess
import shutil

def get_gpu_count():
    """
    Detects the number of available NVIDIA GPUs using nvidia-smi.
    Returns 0 if nvidia-smi is not found or fails.
    """
    if not shutil.which("nvidia-smi"):
        return 0
    
    try:
        # Run nvidia-smi to list GPUs
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=count", "--format=csv,noheader"],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            output = result.stdout.strip()
            if output:
                # The output is usually the count repeated for each GPU, or just the count.
                # Actually --query-gpu=count returns the count *per line*? No, it returns the count.
                # Let's try listing UUIDs and counting lines, it's safer.
                pass
        
        # Safer method: List UUIDs
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=uuid", "--format=csv,noheader"],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            return len([l for l in lines if l.strip()])
            
        return 0
    except Exception:
        return 0

import os
import time
import signal
import atexit

class OllamaCluster:
    def __init__(self, base_port=11500):
        self.base_port = base_port
        self.processes = []
        self.urls = []
        self.gpu_count = get_gpu_count()

    def start(self):
        """
        Starts an Ollama instance for each detected GPU.
        Returns the list of URLs.
        """
        if self.gpu_count <= 1:
            print("Single GPU or no GPU detected. Skipping cluster launch.")
            return []

        print(f"Starting Ollama cluster for {self.gpu_count} GPUs...")
        
        # Locate models
        home = os.path.expanduser("~")
        orig_models = os.path.join(home, ".ollama", "models")
        
        # Check if user models exist
        if not os.path.exists(os.path.join(orig_models, "manifests")) or not os.listdir(os.path.join(orig_models, "manifests")):
            # Check system models
            system_models = "/usr/share/ollama/.ollama/models"
            if os.path.exists(system_models):
                print(f"Using system-wide models from {system_models}")
                orig_models = system_models
            else:
                print(f"Warning: No models found in {orig_models} or {system_models}.")

        for i in range(self.gpu_count):
            port = self.base_port + i
            
            # Setup isolated environment
            iso_home = os.path.join(home, ".ollama_cluster", f"gpu{i}")
            iso_models_dir = os.path.join(iso_home, ".ollama", "models")
            
            os.makedirs(os.path.join(iso_home, ".ollama"), exist_ok=True)
            
            # Symlink models
            if os.path.exists(orig_models):
                if os.path.exists(iso_models_dir):
                    if os.path.islink(iso_models_dir):
                        os.unlink(iso_models_dir)
                    elif os.path.isdir(iso_models_dir):
                        shutil.rmtree(iso_models_dir)
                
                try:
                    os.symlink(orig_models, iso_models_dir)
                except OSError as e:
                    print(f"Warning: Could not symlink models for GPU {i}: {e}")

            # Prepare environment variables
            env = os.environ.copy()
            env["CUDA_VISIBLE_DEVICES"] = str(i)
            env["OLLAMA_HOST"] = f"0.0.0.0:{port}"
            env["HOME"] = iso_home
            
            # Launch process
            print(f"  Launching Ollama on GPU {i} (Port {port})...")
            try:
                # Redirect output to log files
                log_file = open(f"ollama_gpu{i}.log", "w")
                proc = subprocess.Popen(
                    ["ollama", "serve"],
                    env=env,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    preexec_fn=os.setsid # Create new process group for clean kill
                )
                self.processes.append((proc, log_file))
                self.urls.append(f"http://localhost:{port}")
            except Exception as e:
                print(f"  Failed to launch on GPU {i}: {e}")

        # Wait a bit for startup
        print("Waiting for cluster to initialize...")
        time.sleep(2)
        
        return self.urls

    def stop(self):
        """
        Stops all started Ollama processes.
        """
        if not self.processes:
            return

        print("\nStopping Ollama cluster...")
        for proc, log_file in self.processes:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                proc.wait(timeout=5)
            except Exception:
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except:
                    pass
            
            if log_file:
                log_file.close()
        
        self.processes = []
        print("Cluster stopped.")

def get_ollama_launch_commands(gpu_count, base_port=11434):
    """
    Generates a list of shell commands to launch Ollama instances for each GPU.
    """
    commands = []
    for i in range(gpu_count):
        port = base_port + i
        cmd = f"CUDA_VISIBLE_DEVICES={i} OLLAMA_HOST=0.0.0.0:{port} ollama serve"
        commands.append(cmd)
    return commands
