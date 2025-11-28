"""
MAInr - SRA Data Mining Agent

Main entry point for the MAInr application.
"""

import argparse
import atexit
import os
import sys
from pathlib import Path

import pandas as pd
from Bio import Entrez

from src.processing.pipeline import Pipeline
from src.utils.logger import setup_logger
from src.utils.gpu import get_gpu_count, OllamaCluster
from src.utils.hardware import get_gpu_info, estimate_max_workers
from config.settings import ENTREZ_EMAIL, ENTREZ_API_KEY
from src.constants import (
    ValidationLimits,
    PipelineDefaults,
    FileDefaults
)

# Setup main logger
logger = setup_logger("MAInr")


def validate_workers(num_workers: int) -> int:
    """
    Validate and clamp number of workers to acceptable range.

    Args:
        num_workers: Requested number of workers

    Returns:
        Validated number of workers
    """
    if num_workers < ValidationLimits.MIN_WORKERS:
        logger.warning(
            f"Workers {num_workers} below minimum, using {ValidationLimits.MIN_WORKERS}"
        )
        return ValidationLimits.MIN_WORKERS

    if num_workers > ValidationLimits.MAX_WORKERS:
        logger.warning(
            f"Workers {num_workers} above maximum, capping at {ValidationLimits.MAX_WORKERS}"
        )
        return ValidationLimits.MAX_WORKERS

    return num_workers


def validate_ollama_threads(num_threads: int) -> int:
    """
    Validate Ollama thread count.

    Args:
        num_threads: Requested thread count

    Returns:
        Validated thread count
    """
    if num_threads < ValidationLimits.MIN_OLLAMA_THREADS:
        logger.warning(
            f"Ollama threads {num_threads} below minimum, using {ValidationLimits.MIN_OLLAMA_THREADS}"
        )
        return ValidationLimits.MIN_OLLAMA_THREADS

    if num_threads > ValidationLimits.MAX_OLLAMA_THREADS:
        logger.warning(
            f"Ollama threads {num_threads} above maximum, capping at {ValidationLimits.MAX_OLLAMA_THREADS}"
        )
        return ValidationLimits.MAX_OLLAMA_THREADS

    return num_threads


def main():
    """Main application entry point."""
    parser = argparse.ArgumentParser(
        description="MAInr - SRA Mining Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s "drought stress in tomato"
  %(prog)s "heat stress in arabidopsis" -O ./results -n 20
  %(prog)s "cold stress in maize" -t 4 --target-projects 500
        """
    )

    parser.add_argument(
        "topic",
        nargs="?",
        help="Research topic (e.g., 'drought stress in tomato')"
    )
    parser.add_argument(
        "-O", "--output-dir",
        default=FileDefaults.DEFAULT_OUTPUT_DIR,
        help=f"Output directory for results (default: {FileDefaults.DEFAULT_OUTPUT_DIR})"
    )
    parser.add_argument(
        "-n", "--num-workers",
        type=int,
        default=PipelineDefaults.DEFAULT_NUM_WORKERS,
        help=f"Number of threads for parallel processing (default: {PipelineDefaults.DEFAULT_NUM_WORKERS})"
    )
    parser.add_argument(
        "-t", "--ollama-threads",
        type=int,
        default=None,
        help="Number of threads for Ollama (optional)"
    )
    parser.add_argument(
        "--email",
        default=ENTREZ_EMAIL,
        help="Email for NCBI Entrez (or set ENTREZ_EMAIL env var)"
    )
    parser.add_argument(
        "--api-key",
        default=ENTREZ_API_KEY,
        help="API Key for NCBI Entrez (or set ENTREZ_API_KEY env var)"
    )
    parser.add_argument(
        "--target-projects",
        type=int,
        default=None,
        help="Target number of unique projects to retrieve (optional)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        import logging
        logger.setLevel(logging.DEBUG)
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)

    logger.info("=" * 60)
    logger.info("MAInr - SRA Mining Agent")
    logger.info("=" * 60)

    # Configure Entrez
    if args.email:
        Entrez.email = args.email
        logger.debug(f"NCBI Entrez email configured: {args.email}")
    else:
        logger.error("Email is required to use NCBI Entrez")
        logger.error("  Use the --email argument or set the ENTREZ_EMAIL environment variable")
        return 1

    if args.api_key:
        Entrez.api_key = args.api_key
        logger.debug("NCBI Entrez API key configured")

    # Get research topic
    topic = args.topic
    if not topic:
        try:
            topic = input("Enter research topic (e.g., 'drought stress in tomato'): ")
            if not topic:
                topic = "drought stress in tomato"
                logger.info(f"Using default topic: {topic}")
        except (KeyboardInterrupt, EOFError):
            logger.info("\nOperation cancelled by user")
            return 0

    logger.info(f"Research topic: {topic}")

    # Create output directory
    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created output directory: {output_dir}")
        except Exception as e:
            logger.error(f"Failed to create output directory: {e}")
            return 1

    # --- Multi-GPU Cluster Automation ---
    cluster = None
    if not os.getenv("OLLAMA_URL"):
        gpu_count = get_gpu_count()
        if gpu_count > 1:
            logger.info(f"Detected {gpu_count} GPUs - launching optimized Ollama cluster")
            cluster = OllamaCluster()
            urls = cluster.start()

            if urls:
                # Set environment variable for this process
                os.environ["OLLAMA_URL"] = ",".join(urls)
                logger.info(f"Cluster active with {len(urls)} instance(s)")

                # Register cleanup
                atexit.register(cluster.stop)
            else:
                logger.warning("Failed to start cluster, using default configuration")

    # --- Hardware Detection and Configuration ---
    logger.info("")
    logger.info("--- Hardware Detection ---")

    gpus = get_gpu_info()
    num_workers = args.num_workers

    if gpus:
        total_free_mem = sum(gpu['free_memory_mb'] for gpu in gpus)

        logger.info(f"Detected {len(gpus)} GPU(s):")
        for i, gpu in enumerate(gpus):
            logger.info(
                f"  GPU {i}: {gpu['name']} | "
                f"Free VRAM: {gpu['free_memory_mb']} MB / {gpu['total_memory_mb']} MB"
            )

        # Estimate optimal workers
        recommended_workers = estimate_max_workers(gpus, model_size_gb=2.5)

        # Ask user for confirmation
        try:
            user_input = input(
                f"\nEnter number of workers to use [Default: {recommended_workers}]: "
            )
            if user_input.strip():
                try:
                    num_workers = int(user_input)
                except ValueError:
                    logger.warning("Invalid input, using recommended value")
                    num_workers = recommended_workers
            else:
                num_workers = recommended_workers
        except (KeyboardInterrupt, EOFError):
            logger.info("\nUsing recommended value")
            num_workers = recommended_workers
    else:
        logger.info("No NVIDIA GPUs detected (CPU mode)")

    # Validate inputs
    num_workers = validate_workers(num_workers)
    logger.info(f"Using {num_workers} concurrent workers")

    if args.ollama_threads:
        ollama_threads = validate_ollama_threads(args.ollama_threads)
        logger.info(f"Ollama threads: {ollama_threads}")
    else:
        ollama_threads = None

    logger.info(
        "Tip: Set OLLAMA_NUM_PARALLEL to at least this number "
        "for maximum throughput in your Ollama server"
    )

    # --- Run Pipeline ---
    logger.info("")
    logger.info("--- Starting Pipeline ---")

    try:
        pipeline = Pipeline(
            ollama_threads=ollama_threads,
            target_projects=args.target_projects
        )
        results = pipeline.run(topic, max_workers=num_workers)

        # --- Save Results ---
        if results:
            logger.info("")
            logger.info("--- Saving Results ---")

            df = pd.DataFrame(results)
            filename = f"MAInr_results_{topic.replace(' ', '_')}.csv"
            output_path = output_dir / filename

            df.to_csv(output_path, index=False)

            logger.info(f"Results saved to: {output_path}")
            logger.info(f"Total projects analyzed: {len(results)}")

            # Show preview
            if len(df) > 0:
                logger.info("\nPreview of results:")
                preview_cols = ['bioproject', 'title', 'summary']
                available_cols = [col for col in preview_cols if col in df.columns]
                print(df[available_cols].head())

            return 0
        else:
            logger.warning("No results generated")
            return 1

    except KeyboardInterrupt:
        logger.info("\n\nOperation cancelled by user")
        return 130  # Standard exit code for SIGINT
    except Exception as e:
        logger.error(f"Pipeline failed with error: {type(e).__name__}: {e}", exc_info=True)
        return 1
    finally:
        # Cleanup is handled by atexit
        pass


if __name__ == "__main__":
    sys.exit(main() or 0)
