import argparse
import os
import pandas as pd
from Bio import Entrez
from src.processing.pipeline import Pipeline
from config.settings import ENTREZ_EMAIL, ENTREZ_API_KEY

def main():
    parser = argparse.ArgumentParser(description="MAInr - SRA Mining Agent")
    parser.add_argument("topic", nargs="?", help="Research topic (e.g., 'drought stress in tomato')")
    parser.add_argument("-O", "--output-dir", default=".", help="Output directory for results")
    parser.add_argument("-n", "--num-workers", type=int, default=15, help="Number of threads for parallel processing (default: 15)")
    parser.add_argument("-t", "--ollama-threads", type=int, default=None, help="Number of threads for Ollama (optional)")
    parser.add_argument("--email", default=ENTREZ_EMAIL, help="Email for NCBI Entrez")
    parser.add_argument("--api-key", default=ENTREZ_API_KEY, help="API Key for NCBI Entrez")
    
    args = parser.parse_args()

    print("Starting MAInr - SRA Mining Agent")
    
    # Configure Entrez
    if args.email:
        Entrez.email = args.email
    
    if args.api_key:
        Entrez.api_key = args.api_key
        
    if not Entrez.email:
        print("Error: An email is required to use NCBI Entrez.")
        print("  Use the --email argument or set the ENTREZ_EMAIL environment variable.")
        return
    
    topic = args.topic
    if not topic:
        topic = input("Enter research topic (e.g., 'drought stress in tomato'): ")
        if not topic:
            topic = "drought stress in tomato" # Default

    # Create output directory if it doesn't exist
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"Directory created: {args.output_dir}")
    
    # Hardware Detection and Configuration
    from src.utils.hardware import get_gpu_info, estimate_max_workers
    
    print("\n--- Hardware Detection ---")
    gpus = get_gpu_info()
    num_workers = args.num_workers
    
    if gpus:
        print(f"Detected {len(gpus)} GPUs:")
        total_free_mem = 0
        for i, gpu in enumerate(gpus):
            print(f"  GPU {i}: {gpu['name']} | Free VRAM: {gpu['free_memory_mb']} MB / {gpu['total_memory_mb']} MB")
            total_free_mem += gpu['free_memory_mb']
            
        # Estimate based on ~2GB per concurrent stream (assuming model is loaded and shared, mainly KV cache cost)
        # This is a heuristic. Qwen2.5:14b needs ~9GB base.
        # If we assume model is loaded once, extra streams are cheaper.
        recommended_workers = estimate_max_workers(gpus, model_size_gb=2.5) 
        
        print(f"\nTotal Free VRAM: {total_free_mem} MB")
        print(f"Recommended concurrent workers (threads): {recommended_workers}")
        
        try:
            user_input = input(f"Enter number of workers to use [Default: {recommended_workers}]: ")
            if user_input.strip():
                num_workers = int(user_input)
            else:
                num_workers = recommended_workers
        except ValueError:
            print("Invalid input. Using default.")
            num_workers = recommended_workers
    else:
        print("No NVIDIA GPUs detected via nvidia-smi. Using CPU or default settings.")
        
    print(f"Using {num_workers} concurrent workers.")

    pipeline = Pipeline(ollama_threads=args.ollama_threads)
    results = pipeline.run(topic, max_workers=num_workers)
    
    if results:
        df = pd.DataFrame(results)
        filename = f"MAInr_results_{topic.replace(' ', '_')}.csv"
        output_path = os.path.join(args.output_dir, filename)
        df.to_csv(output_path, index=False)
        print(f"\nDone! Results saved to {output_path}")
        print(df[['bioproject', 'title', 'summary']].head())
    else:
        print("\nNo results generated.")

if __name__ == "__main__":
    main()
