import argparse
import os
import pandas as pd
from Bio import Entrez
from src.processing.pipeline import Pipeline
from config.settings import ENTREZ_EMAIL, ENTREZ_API_KEY

def main():
    parser = argparse.ArgumentParser(description="MAInr - Agente de Minería SRA")
    parser.add_argument("topic", nargs="?", help="Tema de investigación (ej., 'drought stress in tomato')")
    parser.add_argument("-O", "--output-dir", default=".", help="Directorio de salida para los resultados")
    parser.add_argument("-n", "--num-workers", type=int, default=15, help="Número de hilos para procesamiento paralelo (default: 15)")
    parser.add_argument("-t", "--ollama-threads", type=int, default=None, help="Número de hilos para Ollama (opcional)")
    parser.add_argument("--email", default=ENTREZ_EMAIL, help="Correo electrónico para NCBI Entrez")
    parser.add_argument("--api-key", default=ENTREZ_API_KEY, help="Clave API para NCBI Entrez")
    
    args = parser.parse_args()

    print("Iniciando MAInr - Agente de Mineria SRA")
    
    # Configurar Entrez
    if args.email:
        Entrez.email = args.email
    
    if args.api_key:
        Entrez.api_key = args.api_key
        
    if not Entrez.email:
        print("Error: Se requiere un correo electronico para usar NCBI Entrez.")
        print("  Usa el argumento --email o establece la variable de entorno ENTREZ_EMAIL.")
        return
    
    topic = args.topic
    if not topic:
        topic = input("Introduce el tema de investigación (ej., 'drought stress in tomato'): ")
        if not topic:
            topic = "drought stress in tomato" # Por defecto

    # Crear directorio de salida si no existe
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print("Directorio creado: {args.output_dir}")
    
    pipeline = Pipeline(ollama_threads=args.ollama_threads)
    results = pipeline.run(topic, max_workers=args.num_workers)
    
    if results:
        df = pd.DataFrame(results)
        filename = f"MAInr_results_{topic.replace(' ', '_')}.csv"
        output_path = os.path.join(args.output_dir, filename)
        df.to_csv(output_path, index=False)
        print(f"\nListo! Resultados guardados en {output_path}")
        print(df[['bioproject', 'title', 'summary']].head())
    else:
        print("\nNo se generaron resultados.")

if __name__ == "__main__":
    main()
