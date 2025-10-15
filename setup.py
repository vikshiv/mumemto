import os
import subprocess
import sys
from setuptools import setup, find_packages, Command
from setuptools.command.build_py import build_py
from setuptools.command.install import install
import shutil
import site

class CMakeBuild(Command):
    description = 'Build C++ extension using CMake'
    user_options = []

    def initialize_options(self):
        self.build_temp = 'build'
        self.build_lib = None
        self.install_scripts = None

    def finalize_options(self):
        self.build_lib = self.get_finalized_command('build').build_lib
        self.install_scripts = self.get_finalized_command('install').install_scripts

    def run(self):
        # Create build directory
        os.makedirs(self.build_temp, exist_ok=True)
        
        subprocess.check_call(['cmake', '-S', '.', '-B', self.build_temp, '-DCMAKE_POLICY_VERSION_MINIMUM=3.5'])
        subprocess.check_call(['cmake', '--build', self.build_temp])
        subprocess.check_call(['cmake', '--install', self.build_temp])
        
        # Get the package directory and scripts directory
        pkg_dir = os.path.join(self.build_lib, 'mumemto')
        bin_dir = self.install_scripts
        
        os.makedirs(pkg_dir, exist_ok=True)
        os.makedirs(bin_dir, exist_ok=True)
        
        # Copy executables to scripts directory and make them executable
        executables = [
            (os.path.join(self.build_temp, 'mumemto_exec'), 'mumemto_exec'),
            (os.path.join(self.build_temp, 'compute_lengths'), 'compute_lengths'),
            (os.path.join(self.build_temp, 'extract_mums'), 'extract_mums'),
            (os.path.join(self.build_temp, 'anchor_merge'), 'anchor_merge'),
        ]
        
        for src, dst in executables:
            dst_path = os.path.join(bin_dir, dst)
            shutil.copy2(src, dst_path)
            os.chmod(dst_path, 0o755)  # Make executable
        
        # Copy the Python script separately
        shutil.copy2(
            'mumemto/mumemto',
            os.path.join(bin_dir, 'mumemto')
        )
        os.chmod(os.path.join(bin_dir, 'mumemto'), 0o755)
        
        # Copy Python files
        for py_file in os.listdir('mumemto'):
            if py_file.endswith('.py'):
                shutil.copy2(
                    os.path.join('mumemto', py_file),
                    os.path.join(pkg_dir, py_file)
                )
        
        # Create __init__.py
        open(os.path.join(pkg_dir, '__init__.py'), 'a').close()

class CustomBuildPy(build_py):
    def run(self):
        self.run_command('cmake_build')
        super().run()

class CustomInstall(install):
    def run(self):
        self.run_command('cmake_build')
        super().run()

def read_requirements():
    requirements = [
        'matplotlib',
        'numpy',
        'tqdm',
        'numba'
    ]
    return requirements

setup(
    name="mumemto",
    version="1.3.2",
    packages=find_packages(),
    install_requires=read_requirements(),
    scripts=['mumemto/mumemto'],  # Only include the Python script here
    package_data={
        'mumemto': ['*.py'],
    },
    cmdclass={
        'cmake_build': CMakeBuild,
        'build_py': CustomBuildPy,
        'install': CustomInstall,
    },
    description="Finding maximal unique matches across pangenomes",
    long_description="Mumemto is a tool for finding a variety of matches across collections of sequences like a pangenome. It includes a visualization tool for visualizing pangenome synteny.",
    author="vikshiv",
    url="https://github.com/vikshiv/mumemto",
    license="GPL-3.0-only",
    python_requires='>=3.6',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
) 
