# MAInr - SRA Data Mining Agent

MAInr (Mining Agent for SRA) is an advanced tool designed to search, extract, and analyze metadata from the NCBI Sequence Read Archive (SRA) database. It uses Large Language Models (LLMs) to generate intelligent search strategies and analyze the context of found studies.

## Features

*   **Intelligent Search**: Generates complex, high-recall search queries using LLM models (Mistral/Qwen) to ensure relevant studies are not missed.
*   **Deep Metadata Extraction**: Retrieves full study details including abstracts and experimental designs using direct XML parsing, providing rich context for analysis.
*   **Massive Mining**: Capable of processing thousands of SRA records, handling pagination and data retrieval efficiently.
*   **Contextual Analysis**: Analyzes metadata for each BioProject using LLMs to extract key information such as experimental conditions, studied tissues, study design, and relevance.
*   **Multi-GPU Cluster Automation**: Automatically detects multiple GPUs and launches an optimized, load-balanced Ollama cluster (one instance per GPU) for maximum throughput.
*   **Parallel Processing**: Uses threads to analyze multiple projects simultaneously, speeding up the mining process.
*   **Hardware Detection**: Automatically detects available GPUs and recommends optimal concurrency settings.
*   **Structured Output**: Generates CSV files with all collected and analyzed information ready for use.
*   **Professional Logging**: Comprehensive logging system with colored console output and file logging for debugging.
*   **Progress Tracking**: Real-time progress bars using tqdm for better visibility into pipeline status.
*   **Robust Error Handling**: Comprehensive exception handling with automatic retries and graceful degradation.
*   **Type Safety**: Full type hints throughout the codebase for better IDE support and fewer runtime errors.

## Requirements

*   Python 3.8+
*   Internet access (to connect to NCBI Entrez)
*   An Ollama instance running locally (or accessible via network) with the configured model (default `qwen2.5:14b`).

## Installation

1.  Clone this repository:
    ```bash
    git clone <repository-url>
    cd MAInr
    ```

2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

3.  Configure Ollama:
    Make sure Ollama is installed and running.
    ```bash
    ollama serve
    ```
    And download the necessary model:
    ```bash
    ollama pull qwen2.5:14b
    ```

## Configuration

### NCBI Entrez Credentials

To use the NCBI API, you need to provide an email (required) and optionally an API Key. You can do this in two ways:

1.  **Environment Variables**:
    Set the `ENTREZ_EMAIL` and `ENTREZ_API_KEY` variables in your system.
    ```bash
    export ENTREZ_EMAIL="your_email@example.com"
    export ENTREZ_API_KEY="your_api_key"
    ```

2.  **Command Line Arguments**:
    Pass them when running the script:
    ```bash
    python3 mainr.py --email "your_email@example.com" --api-key "your_api_key" ...
    ```

### Other Configurations

The `config/settings.py` file contains other default configurations that can be overridden by environment variables or by modifying the file:

*   **OLLAMA_MODEL**: The Ollama model to use (default: `qwen2.5:14b`).
*   **OLLAMA_URL**: URL of the Ollama instance (default: `http://localhost:11434`).
*   **TARGET_PROJECTS**: Target number of unique projects to retrieve.
*   **BATCH_SIZE**: Batch size for SRA searches.

## Usage

Execute the main script from the terminal passing the necessary arguments:

```bash
python3 mainr.py "research topic" [options]
```

### Arguments

*   `topic`: The research topic (e.g., "drought stress in tomato"). If not provided, the system will ask interactively.
*   `-O`, `--output-dir`: Directory where results will be saved (default: current directory).
*   `-n`, `--num-workers`: Number of threads for parallel analysis processing (default: 15).
*   `-t`, `--ollama-threads`: Number of threads Ollama will use for inference (optional).

### Examples

**Basic usage:**
```bash
python3 mainr.py "drought stress in tomato"
```

**Specifying output directory and more worker threads:**
```bash
python3 mainr.py "heat stress in arabidopsis" -O ./results -n 20
```

**Controlling Ollama threads:**
```bash
python3 mainr.py "cold stress in maize" -t 4
```

## Project Structure

*   `mainr.py`: Application entry point (CLI).
*   `config/`: Configuration files.
*   `src/`: Source code.
    *   `llm/`: Client for interacting with Ollama.
    *   `processing/`: Main processing pipeline logic.
    *   `sra/`: Functions to search and retrieve data from NCBI SRA.
    *   `utils/`: Various utilities, including the XML parser.

## Output Columns

The generated CSV file contains the following key columns:

*   **bioproject**: The BioProject accession number (e.g., PRJNA123456).
*   **title**: Title of the study.
*   **summary**: A concise 1-sentence summary of the study's objective (generated by LLM).
*   **organism**: Scientific name of the organism.
*   **study_design**: Classification of the study design (e.g., "Case-Control", "Time Series").
*   **stress_condition**: Specific stress or condition studied.
*   **sequencing_target**: Molecule sequenced (e.g., "mRNA", "Small RNA").
*   **relevance_score**: Relevance score (1-10) for the topic.
*   **tissues**: List of plant tissues identified in the study (e.g., `['root', 'leaf']`).
*   **cultivars**: List of cultivars or genotypes identified (e.g., `['M82']`).
*   **time_course_days**: Maximum duration of the study in days (0 if not a time series).
*   **sra_experiment_count**: Total number of SRA Experiments associated with the BioProject.
*   **biosample_count**: Total number of BioSamples associated with the BioProject.
*   **library_strategy**: Sequencing strategy (e.g., RNA-Seq).
*   **platform**: Sequencing platform (e.g., ILLUMINA).
*   **total_spots**: Total number of reads/spots.
*   **total_bases**: Total number of bases sequenced.

You can find a sample output file here: [examples/example_output.csv](examples/example_output.csv)

## Recent Improvements

This codebase has undergone comprehensive improvements focused on code quality, reliability, and user experience:

### Code Quality
- ✅ Professional logging system with colored output
- ✅ Full type hints across all modules
- ✅ Comprehensive docstrings (Google style)
- ✅ Centralized constants module
- ✅ Improved error handling with specific exception types
- ✅ Progress bars for long-running operations

### Reliability
- ✅ Input validation and safe value clamping
- ✅ Health checks for Ollama cluster
- ✅ Automatic retry logic with exponential backoff
- ✅ Better configuration with reasonable defaults

### Developer Experience
- ✅ Versioned dependencies for reproducibility
- ✅ Additional test suites (LLM client, hardware detection)
- ✅ Comprehensive documentation (CHANGELOG.md, IMPROVEMENTS.md)
- ✅ Better CLI with examples and verbose mode

See [IMPROVEMENTS.md](IMPROVEMENTS.md) and [CHANGELOG.md](CHANGELOG.md) for detailed information.

## Troubleshooting

### Common Issues

**Issue: "Email is required to use NCBI Entrez"**
- Solution: Set the `ENTREZ_EMAIL` environment variable or use `--email` flag

**Issue: "All LLM hosts failed"**
- Solution: Ensure Ollama is running (`ollama serve`) and the model is pulled (`ollama pull qwen2.5:14b`)
- Check logs in `logs/mainr.log` for detailed error messages

**Issue: "No GPUs detected"**
- Solution: This is normal if running on CPU. The tool will work but may be slower.
- If you have GPUs, ensure `nvidia-smi` is available

**Verbose Mode:**
```bash
python3 mainr.py "your topic" -v
```

Check `logs/mainr.log` for detailed execution logs.

## Contribution

Contributions are welcome! Please open an issue or send a pull request for improvements or corrections.

## License

[MIT](LICENSE)
