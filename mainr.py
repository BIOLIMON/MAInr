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
    
    pipeline = Pipeline(ollama_threads=args.ollama_threads)
    results = pipeline.run(topic, max_workers=args.num_workers)
    
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
