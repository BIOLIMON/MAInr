#!/usr/bin/env python3
import pandas as pd
from src.processing.pipeline import Pipeline

def main():
    print("🚀 Iniciando MAInr - Agente de Minería SRA")
    
    # En una aplicación real, esto podría venir de argumentos CLI o un archivo de configuración
    # Por ahora, preguntaremos al usuario o usaremos un valor por defecto para pruebas
    topic = input("Introduce el tema de investigación (ej., 'drought stress in tomato'): ")
    if not topic:
        topic = "drought stress in tomato" # Por defecto
    
    pipeline = Pipeline()
    results = pipeline.run(topic)
    
    if results:
        df = pd.DataFrame(results)
        filename = f"MAInr_results_{topic.replace(' ', '_')}.csv"
        df.to_csv(filename, index=False)
        print(f"\n🎉 ¡Listo! Resultados guardados en {filename}")
        print(df[['bioproject', 'title', 'summary']].head())
    else:
        print("\n⚠️ No se generaron resultados.")

if __name__ == "__main__":
    main()
