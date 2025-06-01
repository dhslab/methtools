import logging
import subprocess
from pathlib import Path
from typing import List, Union, Optional
import tempfile
import os
import shutil # For cleanup of temporary directories

# Third-party import for version parsing.
# Install with: pip install packaging
try:
    from packaging.version import parse as parse_version, InvalidVersion
except ImportError:
    # Provide a fallback or a clear error message if 'packaging' is not installed.
    print("Error: 'packaging' library not found. Please install it using 'pip install packaging'")
    exit(1)

# Configure a basic logger for standalone script testing
# In the actual application, this would be handled by a call to log_config.setup_logging()
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(module)s - %(message)s")

# Get a logger for this module
logger = logging.getLogger(__name__)

def check_tool_version(tool_name: str, min_version_str: str, version_command_flag: str = "--version") -> None:
    """
    Checks if a command-line tool is installed and meets a minimum version.

    Args:
        tool_name: The name of the command-line tool (e.g., "samtools").
        min_version_str: The minimum required version string (e.g., "1.21").
        version_command_flag: The argument to get the tool's version (e.g., "--version").

    Raises:
        RuntimeError: If the tool is not found, its version can't be parsed,
                      or it doesn't meet the minimum version requirement.
    """
    logger.info(f"Checking version for '{tool_name}' (minimum required: {min_version_str})...")
    try:
        min_version = parse_version(min_version_str)
    except InvalidVersion:
        logger.error(f"Invalid minimum version string format provided: '{min_version_str}'")
        raise RuntimeError(f"Invalid minimum version string format: '{min_version_str}'")

    try:
        result = subprocess.run(
            [tool_name, version_command_flag],
            capture_output=True,
            text=True,
            check=True,  # Raise an exception for non-zero exit codes
            encoding="utf-8"
        )
        
        output_lines = result.stdout.strip().splitlines()
        version_str_extracted = None
        for line in output_lines:
            for part in line.split():
                try:
                    current_version_parsed = parse_version(part)
                    if current_version_parsed.major is not None and current_version_parsed.minor is not None:
                        version_str_extracted = str(current_version_parsed)
                        break
                except InvalidVersion:
                    continue
            if version_str_extracted:
                break
        
        if not version_str_extracted:
            logger.error(f"Could not extract a valid version string for '{tool_name}' from output:\n{result.stdout.strip()}")
            raise RuntimeError(f"Could not determine version for '{tool_name}'.")

        current_version = parse_version(version_str_extracted)
        logger.info(f"Found {tool_name} version: {current_version}")

        if current_version < min_version:
            msg = f"'{tool_name}' version {current_version} is older than the required {min_version}."
            logger.error(msg)
            raise RuntimeError(msg)
        
        logger.info(f"'{tool_name}' version {current_version} meets the requirement (>= {min_version}). ✅")

    except FileNotFoundError:
        logger.error(f"Tool not found: '{tool_name}' is not installed or not in the system's PATH.")
        raise RuntimeError(f"Tool not found: '{tool_name}'.")
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.strip() if e.stderr else e.stdout.strip()
        logger.error(f"Error checking '{tool_name}' version. Command exited with an error:\n{stderr}")
        raise RuntimeError(f"Failed to get '{tool_name}' version.")
    except Exception as e:
        logger.error(f"An unexpected error occurred while checking '{tool_name}': {e}", exc_info=True)
        raise RuntimeError(f"An unexpected error checking '{tool_name}'.")

def run_command(cmd: Union[str, List[str]], dry_run: bool = False, workdir: Optional[Union[str, Path]] = None) -> None:
    """
    Executes a shell command with logging and error handling.

    Args:
        cmd: The command to execute, either as a single string (with shell=True)
             or a list of arguments (with shell=False).
        dry_run: If True, the command is logged but not executed.
        workdir: The working directory for the command. Defaults to the current directory.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
        RuntimeError: For other execution-related errors.
    """
    shell = isinstance(cmd, str)
    cmd_str = cmd if shell else " ".join(map(str, cmd))
    
    if workdir:
        logger.info(f"Running command in '{workdir}': {cmd_str}")
    else:
        logger.info(f"Running command: {cmd_str}")

    if dry_run:
        logger.info("--> [Dry Run] Command not executed.")
        return

    try:
        process = subprocess.run(
            cmd,
            cwd=str(workdir) if workdir else None,
            shell=shell,
            check=True,
            text=True,
            capture_output=True
        )
        if process.stdout:
            logger.debug(f"Command stdout:\n{process.stdout.strip()}")
        logger.info("--> Command completed successfully. ✅")

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}: {cmd_str}")
        if e.stdout:
            logger.error(f"STDOUT from failed command:\n{e.stdout.strip()}")
        if e.stderr:
            logger.error(f"STDERR from failed command:\n{e.stderr.strip()}")
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred trying to run command '{cmd_str}': {e}", exc_info=True)
        raise RuntimeError(f"Failed to run command '{cmd_str}'.")

def get_temporary_directory(
    user_specified_dir: Optional[Union[str, Path]] = None,
    prefix: str = "methtools_tmp_"
) -> Path:
    """
    Determines, creates, and returns a temporary directory path.

    This function prioritizes directory locations in the following order:
    1. A user-specified directory (if provided).
    2. Common HPC environment variables (e.g., METHTOOLS_TMPDIR, SLURM_TMPDIR, TMPDIR).
    3. The system's default temporary location.
    4. As a last resort, a subdirectory in the current working directory.

    The created directory will have a unique name starting with the given prefix.
    It is the responsibility of the caller to clean up this directory when done.

    Args:
        user_specified_dir: An optional path (str or Path) to a base directory
                            where the temporary directory should be created.
        prefix: A prefix for the temporary directory name.

    Returns:
        A pathlib.Path object pointing to the created temporary directory.

    Raises:
        OSError: If a temporary directory cannot be created in any of the
                 attempted locations.
    """
    base_dirs_to_try: List[Path] = []

    if user_specified_dir:
        user_path = Path(user_specified_dir)
        if not user_path.is_dir():
            try:
                logger.info(f"User-specified temp base directory '{user_path}' does not exist. Attempting to create it.")
                user_path.mkdir(parents=True, exist_ok=True)
                base_dirs_to_try.append(user_path)
            except OSError as e:
                logger.warning(f"Could not create user-specified temp base directory '{user_path}': {e}. Ignoring.")
        else:
            base_dirs_to_try.append(user_path)


    # Check common HPC and general environment variables for temporary directories
    env_vars_priority = [
        "METHTOOLS_TMPDIR",  # Custom override for this specific tool
        "SLURM_TMPDIR",    # For Slurm HPC environments
        "LSF_TMPDIR",      # For LSF HPC environments
        "TMPDIR",          # Common on Unix-like systems
        "TEMP",            # Common on Windows
        "TMP"              # Common on Windows
    ]
    for env_var in env_vars_priority:
        env_val = os.environ.get(env_var)
        if env_val:
            env_path = Path(env_val)
            if env_path.is_dir(): # Only consider if the path from env var actually exists and is a directory
                base_dirs_to_try.append(env_path)
                # Optional: break here to use the first valid env var found
                # break 
    
    # Attempt to create the temporary directory in one of the prioritized base locations
    for base_path in base_dirs_to_try:
        try:
            # tempfile.mkdtemp creates a directory with a unique name
            temp_dir_path_str = tempfile.mkdtemp(prefix=prefix, dir=base_path)
            logger.info(f"Using temporary directory: {temp_dir_path_str} (based on: {base_path})")
            return Path(temp_dir_path_str)
        except OSError as e:
            logger.warning(f"Could not create temp directory in '{base_path}': {e}. Trying next option.")
            continue # Try the next base directory
    
    # If all specified/environment variable locations fail, fallback to system default
    try:
        temp_dir_path_str = tempfile.mkdtemp(prefix=prefix)
        logger.info(f"Using system default temporary directory: {temp_dir_path_str}")
        return Path(temp_dir_path_str)
    except OSError as e_sys:
        logger.warning(f"Failed to create temporary directory in system default location: {e_sys}.")
        # As a last resort, try creating in a subdirectory of the current working directory
        try:
            local_fallback_base = Path("./").resolve() / f".{prefix}local_fallback"
            local_fallback_base.mkdir(parents=True, exist_ok=True) # Ensure the base for local fallback exists
            temp_dir_path_str = tempfile.mkdtemp(prefix=prefix, dir=local_fallback_base)
            logger.info(f"Using local fallback temporary directory: {temp_dir_path_str} (within {local_fallback_base})")
            return Path(temp_dir_path_str)
        except OSError as e_local:
            logger.error(f"Critically failed to create any temporary directory, even in local fallback: {e_local}", exc_info=True)
            raise OSError("Unable to create a temporary directory in any specified or default location.") from e_local


# --- Test Suite ---
if __name__ == "__main__":
    print("-" * 50)
    print("--- Running Tests for system.py ---")
    print("-" * 50)

    # --- Tests for check_tool_version ---
    print("\n[1] Testing 'check_tool_version'...")
    try:
        check_tool_version("python3", "3.6", version_command_flag="-V") # Test 1.1
    except Exception as e: logger.error(f"Test 1.1 FAILED: {e}\n")
    try:
        check_tool_version("python3", "99.0", version_command_flag="-V") # Test 1.2
    except RuntimeError as e: logger.info(f"Test 1.2 PASSED as expected: Caught required error -> {e}\n")
    except Exception as e: logger.error(f"Test 1.2 FAILED with unexpected error: {e}\n")
    try:
        check_tool_version("nonexistenttool12345", "1.0") # Test 1.3
    except RuntimeError as e: logger.info(f"Test 1.3 PASSED as expected: Caught required error -> {e}\n")
    except Exception as e: logger.error(f"Test 1.3 FAILED with unexpected error: {e}\n")

    # --- Tests for run_command ---
    print("\n[2] Testing 'run_command'...")
    try: # Test 2.1
        print("--> Testing successful command 'echo'")
        run_command('echo "Hello from run_command!"')
        logger.info("Test 2.1 PASSED\n")
    except Exception as e: logger.error(f"Test 2.1 FAILED: {e}\n")
    try: # Test 2.2
        print("--> Testing dry_run functionality")
        run_command('echo "This should not appear"', dry_run=True)
        logger.info("Test 2.2 PASSED\n")
    except Exception as e: logger.error(f"Test 2.2 FAILED: {e}\n")
    try: # Test 2.3
        print("--> Testing a failing command 'ls non_existent_dir'")
        run_command("ls a_directory_that_does_not_exist_12345")
    except subprocess.CalledProcessError: logger.info("Test 2.3 PASSED as expected: Caught CalledProcessError\n")
    except Exception as e: logger.error(f"Test 2.3 FAILED with unexpected error: {e}\n")
    try: # Test 2.4
        print("--> Testing command as a list")
        run_command(["echo", "Hello", "from", "a", "list", "command"])
        logger.info("Test 2.4 PASSED\n")
    except Exception as e: logger.error(f"Test 2.4 FAILED: {e}\n")

    # --- Tests for get_temporary_directory ---
    print("\n[3] Testing 'get_temporary_directory'...")
    created_temp_dirs = []
    try: # Test 3.1: System default
        print("--> Testing default temporary directory creation")
        default_temp_dir = get_temporary_directory(prefix="test_default_")
        assert default_temp_dir.exists()
        assert default_temp_dir.is_dir()
        assert default_temp_dir.name.startswith("test_default_")
        created_temp_dirs.append(default_temp_dir)
        logger.info(f"Test 3.1 PASSED: Created {default_temp_dir}\n")
    except Exception as e:
        logger.error(f"Test 3.1 FAILED: {e}\n")

    try: # Test 3.2: User-specified existing directory
        print("--> Testing user-specified temporary directory creation (in existing dir)")
        user_base_dir = Path("./my_custom_temp_base").resolve()
        user_base_dir.mkdir(exist_ok=True) # Ensure base exists
        user_temp_dir = get_temporary_directory(user_specified_dir=user_base_dir, prefix="test_user_")
        assert user_temp_dir.exists()
        assert user_temp_dir.is_dir()
        assert user_temp_dir.name.startswith("test_user_")
        assert user_temp_dir.parent == user_base_dir
        created_temp_dirs.append(user_temp_dir)
        logger.info(f"Test 3.2 PASSED: Created {user_temp_dir}\n")
        # Clean up the base for this test if it was created
        if user_base_dir.name == "my_custom_temp_base": # Be careful with rmtree
             shutil.rmtree(user_base_dir) 
    except Exception as e:
        logger.error(f"Test 3.2 FAILED: {e}\n")

    try: # Test 3.3: User-specified non-existing directory (should be created)
        print("--> Testing user-specified temporary directory creation (in non-existing dir)")
        user_base_dir_ne = Path("./my_ne_custom_temp_base").resolve()
        if user_base_dir_ne.exists(): shutil.rmtree(user_base_dir_ne) # ensure it doesn't exist
        
        user_temp_dir_ne = get_temporary_directory(user_specified_dir=user_base_dir_ne, prefix="test_user_ne_")
        assert user_temp_dir_ne.exists()
        assert user_temp_dir_ne.is_dir()
        assert user_temp_dir_ne.name.startswith("test_user_ne_")
        assert user_temp_dir_ne.parent == user_base_dir_ne
        created_temp_dirs.append(user_temp_dir_ne) # Add to list for cleanup
        logger.info(f"Test 3.3 PASSED: Created {user_temp_dir_ne}\n")
    except Exception as e:
        logger.error(f"Test 3.3 FAILED: {e}\n")
    finally:
        # Clean up all created temporary directories
        for d in created_temp_dirs:
            if d.exists():
                shutil.rmtree(d)
                # logger.info(f"Cleaned up test temp dir: {d}")
        # Clean up specific test base dirs if they still exist due to test failure before cleanup
        test_base_dir_ne = Path("./my_ne_custom_temp_base").resolve()
        if test_base_dir_ne.exists(): shutil.rmtree(test_base_dir_ne)
        local_fallback_base = Path("./").resolve() / f".methtools_tmp_local_fallback" # Match prefix in function
        if local_fallback_base.exists(): shutil.rmtree(local_fallback_base)


    print("-" * 50)
    print("--- Test Run Complete ---")
    print("-" * 50)
