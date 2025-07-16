import subprocess
import sys
from pathlib import Path


# Fortran source files
FORTRAN_SOURCES = ["idahopot.f", "av18.f", "cdbonn.f", "lo450.f", "lo500.f", "lo550.f", "nlo450.f", "nlo500.f", "nlo550.f", "n2lo450.f", "n2lo500.f", "n2lo550.f", "n3lo450new.f", "n3lo500new.f", "n3lo550new.f", "n4lo450.f", "n4lo500.f", "n4lo550.f", "nnlo_opt.f", "nnlo_sat.f", "n3lo.f"]

# f2py wrap interface file
INTERFACE_FILE = "wrapforce.pyf"
WRAPPER_F90 = "wrapforce.f90"

# python module name
MODULE_NAME = "force_module"


def run_command(command, error_message):
    print(f"Executing command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        print("error: 'f2py' command not found.", file=sys.stderr)
        print("   Please ensure you have NumPy installed (pip install numpy)", file=sys.stderr)
        print("   and that 'f2py' is in your system PATH.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        # If compilation fails, print detailed stdout and stderr
        print(f"error: {error_message}", file=sys.stderr)
        print("\n--- f2py STDOUT ---\n", e.stdout, file=sys.stderr)
        print("\n--- f2py STDERR ---\n", e.stderr, file=sys.stderr)
        sys.exit(1)


def clean():
    print("--- Cleaning build artifacts...")
    # Added cleaning for .pyf files
    patterns = ["*.o", "*.mod", INTERFACE_FILE, f"{MODULE_NAME}*.so", f"{MODULE_NAME}*.dylib"]
    for pattern in patterns:
        for file in Path.cwd().glob(pattern):
            try:
                file.unlink()
                print(f"  Deleted: {file.name}")
            except OSError as e:
                print(f"  Error: Unable to delete {file.name}: {e}", file=sys.stderr)


def generate_interface():
    print(f"--- Generating interface file: {INTERFACE_FILE}...")
    command = ["f2py", WRAPPER_F90, "-h", INTERFACE_FILE, f"-m{MODULE_NAME}", "--overwrite-signature"]
    run_command(command, f"Failed to generate interface file {INTERFACE_FILE}.")
    print(f"Successfully generated interface file: {INTERFACE_FILE}")


def compile_module():
    """Compile all source files into a Python module using the generated interface file."""
    print(f"--- Compiling module: {MODULE_NAME}...")
    all_sources = [INTERFACE_FILE, WRAPPER_F90] + FORTRAN_SOURCES
    command = ["f2py", "-c", f"-m{MODULE_NAME}", *all_sources]
    run_command(command, f"Failed to compile module {MODULE_NAME}.")
    print(f"Successfully created shared library: {MODULE_NAME}")


def main():
    """Execute the full clean, generate interface, and build process."""
    if not Path(WRAPPER_F90).exists():
        print(f"Error: Main wrapper file {WRAPPER_F90} not found.", file=sys.stderr)
        print("   Please ensure you run this script in the directory containing the Fortran source files.", file=sys.stderr)
        sys.exit(1)

    clean()
    generate_interface()
    compile_module()
    print("\nBuild completed!")


if __name__ == "__main__":
    main()
