# Changelog

All notable changes to MAInr will be documented in this file.

## [Unreleased] - Code Quality Improvements

### Added

- **Professional logging system** (`src/utils/logger.py`)
  - Colored console output for better readability
  - File logging with detailed formatting
  - Configurable log levels
  - Separate log files for different components

- **Constants module** (`src/constants.py`)
  - Centralized configuration values
  - NCBI API limits and constraints
  - Pipeline defaults
  - Ollama configuration
  - GPU estimation parameters
  - Validation limits

- **Health checks for Ollama cluster** (`src/utils/gpu.py`)
  - HTTP-based readiness checks
  - Configurable timeout and retry logic
  - Graceful degradation on failure

- **Progress bars** using tqdm
  - Project collection progress
  - Analysis progress indication
  - Better user experience

- **Comprehensive type hints** across all modules
  - Full typing support for better IDE integration
  - Runtime type safety
  - Improved code documentation

- **Validation system**
  - Input parameter validation
  - Safe value clamping
  - User-friendly error messages

- **Versioned dependencies** (`requirements.txt`)
  - Pinned versions for reproducibility
  - Clear dependency categories
  - Optional dev dependencies

- **Additional tests**
  - `tests/test_llm_client.py`: LLM client tests
  - `tests/test_hardware.py`: Hardware detection tests
  - Better test coverage

### Changed

- **Improved error handling** throughout the codebase
  - Specific exception types instead of bare `except`
  - Better error messages with context
  - Proper logging of errors with stack traces

- **Configuration system** (`config/settings.py`)
  - More reasonable defaults (100 projects instead of 100,000)
  - Better documentation
  - Uses constants for defaults

- **LLM Client** (`src/llm/client.py`)
  - Removed confusing comments about thread safety
  - Better exception handling with specific types
  - Improved logging
  - Full type hints

- **GPU Cluster** (`src/utils/gpu.py`)
  - Added health checks with HTTP polling
  - Better error handling for symlink failures (Windows compatibility)
  - Improved logging
  - Graceful cleanup on errors

- **Hardware Detection** (`src/utils/hardware.py`)
  - Validation for model size parameter
  - Better VRAM calculation with safety margin
  - Improved logging

- **Main Entry Point** (`mainr.py`)
  - Better command-line interface with examples
  - Improved user prompts
  - Verbose mode support (`-v`)
  - Better error handling and exit codes
  - Proper cleanup on interruption

- **Pipeline** (`src/processing/pipeline.py`)
  - Progress bars for long-running operations
  - Better error handling with context
  - Improved logging
  - Type hints throughout

- **SRA Search** (`src/sra/search.py`)
  - Better logging of retry attempts
  - Type hints
  - Improved documentation

### Fixed

- **Critical bugs:**
  - Duplicate import statement in `pipeline.py` - FIXED
  - Incorrect Ollama parameter `num_thread` → `num_threads` - FIXED
  - Hardcoded email in tests - FIXED (now uses env var)
  - Duplicate "LLM Configuration" comment - FIXED

- **Configuration issues:**
  - Unrealistic TARGET_PROJECTS default (100,000 → 100)
  - Unrealistic BATCH_SIZE default (10,000 → 100)
  - Better validation of all user inputs

- **Code quality:**
  - Removed commented-out code in `llm/client.py`
  - Cleaned up unnecessary comments
  - Consistent formatting

- **Windows compatibility:**
  - Fallback to copying when symlinks fail
  - Better handling of file permissions

### Removed

- Unnecessary commented code blocks
- Redundant imports
- Over-complicated thread synchronization comments

## [Previous Versions]

_(Previous version history would go here)_

---

## Migration Guide

### For existing users:

1. **Update dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Configuration changes:**
   - Default `TARGET_PROJECTS` is now 100 (was 100,000)
   - Use `--target-projects` flag to specify custom values
   - Set `ENTREZ_EMAIL` environment variable instead of hardcoding

3. **New environment variables:**
   ```bash
   export ENTREZ_EMAIL="your@email.com"
   export ENTREZ_API_KEY="your_api_key"  # Optional but recommended
   ```

4. **Logging:**
   - Logs are now written to `logs/mainr.log` by default
   - Use `-v` flag for verbose output
   - Check log files for detailed debugging

5. **Progress indication:**
   - Progress bars now show real-time progress
   - Better visibility into pipeline status

### Breaking Changes

- Minimum Python version: 3.8+ (unchanged)
- Some print statements replaced with logging (output format may differ)
- Default configuration values changed (more conservative)

---

## Future Plans

- [ ] YAML configuration file support
- [ ] Result caching system
- [ ] Rate limiting for NCBI API
- [ ] More comprehensive test suite
- [ ] CI/CD with GitHub Actions
- [ ] API documentation with Sphinx
- [ ] Performance profiling and optimization
