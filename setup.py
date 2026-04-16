import glob
import os
import platform
import shutil
import subprocess
import sys
import sysconfig

import pybind11
import numpy
from setuptools import Command, find_packages, setup
from setuptools.command.bdist_wheel import bdist_wheel as bdist_wheel_cmd
from setuptools.command.build_py import build_py
from setuptools.command.install import install

_REPO_ROOT = os.path.dirname(os.path.abspath(__file__))


def _read_readme() -> str:
    path = os.path.join(_REPO_ROOT, "README.md")
    with open(path, encoding="utf-8") as f:
        return f.read()


class PlatformWheel(bdist_wheel_cmd):
    """Binary extensions are built by CMake; wheels must be platform-specific."""

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False


class CMakeBuild(Command):
    description = "Build native code (CLI + libmumemto + pybind11 extension)"
    user_options = []

    def initialize_options(self):
        self.build_temp = "build"
        self.build_lib = None
        self.install_scripts = None

    def finalize_options(self):
        self.build_lib = self.get_finalized_command("build").build_lib
        self.install_scripts = self.get_finalized_command("install").install_scripts

    def _cmake_prefix_path(self):
        return os.path.abspath(os.path.join(self.build_temp, "install"))

    def _cmake_env(self) -> dict:
        """Prefer the active interpreter over other Pythons on PATH (e.g. conda vs base)."""
        env = os.environ.copy()
        bindir = os.path.dirname(sys.executable)
        env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
        return env

    def _cmake_exe(self) -> str:
        bindir = os.path.dirname(sys.executable)
        cand = os.path.join(bindir, "cmake")
        if os.path.isfile(cand):
            return cand
        import shutil

        path = self._cmake_env().get("PATH", "")
        w = shutil.which("cmake", path=path)
        if w:
            return w
        return "cmake"

    def _run_python_bindings_cmake(self, pkg_dir: str) -> None:
        """Configure/build/install _mumemto_core into the package tree."""
        cmake_prefix = self._cmake_prefix_path()
        bindings_src = os.path.join(_REPO_ROOT, "python_bindings")
        bindings_build = os.path.join(self.build_temp, "pybind_build")
        install_prefix = os.path.dirname(pkg_dir)

        os.makedirs(bindings_build, exist_ok=True)

        py_include = sysconfig.get_path("include")
        cx = self._cmake_exe()
        configure = [
            cx,
            "-S",
            bindings_src,
            "-B",
            bindings_build,
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DCMAKE_PREFIX_PATH={cmake_prefix}",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DPython_ROOT_DIR={sys.prefix}",
            f"-DPython_INCLUDE_DIR={py_include}",
            "-DPython_FIND_REGISTRY=NEVER",
            f"-Dpybind11_DIR={pybind11.get_cmake_dir()}",
            f"-DPython_NumPy_INCLUDE_DIR={numpy.get_include()}",
            "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
        ]
        env = self._cmake_env()
        subprocess.check_call(configure, env=env)
        subprocess.check_call([cx, "--build", bindings_build], env=env)
        subprocess.check_call(
            [cx, "--install", bindings_build, "--prefix", install_prefix],
            env=env,
        )

        lib_dir = os.path.join(cmake_prefix, "lib")
        if os.path.isdir(lib_dir):
            for name in os.listdir(lib_dir):
                if name.startswith("libmumemto") and (
                    name.endswith(".so") or ".so." in name
                ):
                    src = os.path.join(lib_dir, name)
                    if os.path.isfile(src):
                        shutil.copy2(src, pkg_dir)
                elif name.startswith("libmumemto") and name.endswith(".dylib"):
                    src = os.path.join(lib_dir, name)
                    if os.path.isfile(src):
                        shutil.copy2(src, pkg_dir)

        # Wheel / local install: depend only on $ORIGIN next to _mumemto_core (Linux).
        if platform.system() == "Linux":
            for path in glob.glob(os.path.join(pkg_dir, "_mumemto_core*.so")):
                try:
                    subprocess.run(
                        ["patchelf", "--set-rpath", "$ORIGIN", path],
                        check=True,
                        capture_output=True,
                    )
                except (FileNotFoundError, subprocess.CalledProcessError):
                    pass

    def run(self):
        os.makedirs(self.build_temp, exist_ok=True)

        cx = self._cmake_exe()
        env = self._cmake_env()
        subprocess.check_call(
            [
                cx,
                "-S",
                ".",
                "-B",
                self.build_temp,
                "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
            ],
            env=env,
        )
        subprocess.check_call([cx, "--build", self.build_temp], env=env)
        subprocess.check_call([cx, "--install", self.build_temp], env=env)

        pkg_dir = os.path.join(self.build_lib, "mumemto")
        bin_dir = self.install_scripts

        os.makedirs(pkg_dir, exist_ok=True)
        os.makedirs(bin_dir, exist_ok=True)

        executables = [
            (os.path.join(self.build_temp, "mumemto_exec"), "mumemto_exec"),
            (os.path.join(self.build_temp, "compute_lengths"), "compute_lengths"),
            (os.path.join(self.build_temp, "extract_mums"), "extract_mums"),
            (os.path.join(self.build_temp, "anchor_merge"), "anchor_merge"),
        ]

        for src, dst in executables:
            dst_path = os.path.join(bin_dir, dst)
            shutil.copy2(src, dst_path)
            os.chmod(dst_path, 0o755)

        shutil.copy2("mumemto/mumemto", os.path.join(bin_dir, "mumemto"))
        os.chmod(os.path.join(bin_dir, "mumemto"), 0o755)

        for py_file in os.listdir("mumemto"):
            if py_file.endswith(".py"):
                shutil.copy2(
                    os.path.join("mumemto", py_file),
                    os.path.join(pkg_dir, py_file),
                )

        self._run_python_bindings_cmake(pkg_dir)


class CustomBuildPy(build_py):
    def run(self):
        self.run_command("cmake_build")
        super().run()


class CustomInstall(install):
    def run(self):
        self.run_command("cmake_build")
        super().run()


def read_requirements():
    return [
        "matplotlib",
        "numpy",
        "tqdm",
        "numba",
    ]


setup(
    name="mumemto",
    version="1.4.0",
    packages=find_packages(),
    install_requires=read_requirements(),
    scripts=["mumemto/mumemto"],
    package_data={
        "mumemto": [
            "*.py",
            "*.so",
            "*.so.*",
            "*.pyd",
            "*.dylib",
        ],
    },
    cmdclass={
        "cmake_build": CMakeBuild,
        "build_py": CustomBuildPy,
        "install": CustomInstall,
        "bdist_wheel": PlatformWheel,
    },
    description="Finding maximal unique matches across pangenomes",
    long_description=_read_readme(),
    long_description_content_type="text/markdown",
    author="vikshiv",
    url="https://github.com/vikshiv/mumemto",
    license="GPL-3.0-only",
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    include_package_data=True,
)
