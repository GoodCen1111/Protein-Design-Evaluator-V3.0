# -*- coding: utf-8 -*-
"""
ProteinDesignEvaluator v3.0 - Build Script
Cross-platform executable builder for Windows and Linux
"""

import os
import sys
import subprocess
import platform
import shutil

def check_dependencies():
    """Check if required packages are installed"""
    required = ['PyQt5', 'numpy', 'matplotlib']
    missing = []

    for package in required:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)

    if missing:
        print(f"Missing packages: {', '.join(missing)}")
        print("Please install them with: pip install -r requirements.txt")
        return False
    return True

def install_pyinstaller():
    """Install PyInstaller if not available"""
    try:
        import PyInstaller
        print(f"PyInstaller version: {PyInstaller.__version__}")
        return True
    except ImportError:
        print("PyInstaller not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])
        return True

def build_windows():
    """Build executable for Windows"""
    print("Building for Windows...")

    # Clean previous builds
    if os.path.exists('dist'):
        shutil.rmtree('dist')
    if os.path.exists('build'):
        shutil.rmtree('build')

    # Build with PyInstaller
    result = subprocess.run(
        ['pyinstaller', 'ProteinDesignEvaluator.spec', '--noconfirm'],
        capture_output=True,
        text=True
    )

    if result.returncode == 0:
        print("Windows build completed successfully!")
        print("Output: dist/ProteinDesignEvaluator/")
        return True
    else:
        print("Build failed:")
        print(result.stderr)
        return False

def build_linux():
    """Build executable for Linux"""
    print("Building for Linux...")

    if os.path.exists('dist'):
        shutil.rmtree('dist')
    if os.path.exists('build'):
        shutil.rmtree('build')

    # Build with PyInstaller (console mode for Linux)
    result = subprocess.run(
        ['pyinstaller', 'main.py', '--name=ProteinDesignEvaluator',
         '--windowed', '--onefile', '--noconfirm'],
        capture_output=True,
        text=True
    )

    if result.returncode == 0:
        print("Linux build completed successfully!")
        print("Output: dist/ProteinDesignEvaluator")
        return True
    else:
        print("Build failed:")
        print(result.stderr)
        return False

def main():
    print("=" * 60)
    print("  ProteinDesignEvaluator v3.0 - Build Script")
    print("=" * 60)
    print()

    # Check Python version
    if sys.version_info < (3, 8):
        print("Error: Python 3.8 or higher is required")
        sys.exit(1)

    # Check dependencies
    if not check_dependencies():
        print()
        response = input("Do you want to install missing packages? (y/n): ")
        if response.lower() == 'y':
            subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
        else:
            sys.exit(1)

    # Install PyInstaller
    install_pyinstaller()

    # Build for current platform
    system = platform.system()
    print(f"Detected system: {system}")
    print()

    if system == "Windows":
        success = build_windows()
    elif system == "Linux":
        success = build_linux()
    else:
        print(f"Unsupported system: {system}")
        print("Building with default settings...")
        success = build_linux() if system == "Linux" else build_windows()

    if success:
        print()
        print("=" * 60)
        print("  Build completed!")
        print("=" * 60)
    else:
        print()
        print("Build failed. Please check the error messages above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
