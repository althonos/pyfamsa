#!/usr/bin/env python3

import configparser
import functools
import glob
import itertools
import io
import multiprocessing.pool
import os
import platform
import re
import setuptools
import setuptools.extension
import subprocess
import sys
import sysconfig
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.extension import Extension
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

try:
    import semantic_version
except ImportError as err:
    semantic_version = err

# --- Constants -----------------------------------------------------------------

SETUP_FOLDER = os.path.relpath(os.path.realpath(os.path.join(__file__, os.pardir)))
FAMSA_FOLDER = os.path.join(SETUP_FOLDER, "vendor", "FAMSA")

def _detect_target_machine(platform):
    if platform == "win32":
        return "x86"
    return platform.rsplit("-", 1)[-1]

def _detect_target_cpu(platform):
    machine = _detect_target_machine(platform)
    if re.match("^mips", machine):
        return "mips"
    elif re.match("^(aarch64|arm64)$", machine):
        return "aarch64"
    elif re.match("^arm", machine):
        return "arm"
    elif re.match("(x86_64)|(x86)|(AMD64|amd64)|(^i.86$)", machine):
        return "x86"
    elif re.match("^(powerpc|ppc)", machine):
        return "ppc"
    return None

def _detect_target_system(platform):
    if platform.startswith("win"):
        return "windows"
    elif platform.startswith("macosx"):
        return "macos"
    elif platform.startswith("linux"):
        return "linux_or_android"
    elif platform.startswith("freebsd"):
        return "freebsd"
    return None


# --- Utils ------------------------------------------------------------------

_HEADER_PATTERN = re.compile("^@@ -(\d+),?(\d+)? \+(\d+),?(\d+)? @@$")


def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


def _patch_osx_compiler(compiler, machine):
    # On newer OSX, Python has been compiled as a universal binary, so
    # it will attempt to pass universal binary flags when building the
    # extension. This will not work because the code makes use of CPU
    # specific SIMD extensions.
    for tool in ("compiler", "compiler_so", "linker_so"):
        flags = getattr(compiler, tool)
        i = next((i for i in range(1, len(flags)) if flags[i-1] == "-arch" and flags[i] != machine), None)
        if i is not None:
            flags.pop(i)
            flags.pop(i-1)
        if "-fPIC" not in flags:
            flags.append("-fPIC")
        while "-O3" in flags:
            i = flags.index("-O3")
            flags[i] = "-O2"


def _apply_patch(s,patch,revert=False):
    # see https://stackoverflow.com/a/40967337
    s = s.splitlines(keepends=True)
    p = patch.splitlines(keepends=True)
    t = []
    i = 0
    sl = 0
    midx, sign = (1,'+') if not revert else (3,'-')
    while i < len(p) and p[i].startswith(("---","+++")):
        i += 1 # skip header lines

    while i < len(p):
        match = _HEADER_PATTERN.match(p[i])
        if not match:
          raise ValueError("Invalid line in patch: {!r}".format(p[i]))
        i += 1
        l = int(match.group(midx)) - 1 + (match.group(midx+1) == '0')
        t.extend(s[sl:l])
        sl = l
        while i < len(p) and p[i][0] != '@':
            if i+1 < len(p) and p[i+1][0] == '\\':
                line = p[i][:-1]
                i += 2
            else:
                line = p[i]
                i += 1
            if len(line) > 0:
                if line[0] == sign or line[0] == ' ':
                    t += line[1:]
                sl += (line[0] != sign)

    t.extend(s[sl:])
    return ''.join(t)


def _is_compiler_clang(compiler):
    try:
        if compiler.compiler_type == "unix":
            proc = subprocess.run([compiler.compiler[0], "--version"], capture_output=True)
            if proc.returncode == 0:
                return any(
                    line.startswith(b"Apple clang")
                    for line in proc.stdout.splitlines()
                )
    except:
        pass
    return False


def _detect_compiler_version(compiler):
    try:
        if compiler.compiler_type == "unix":
            proc = subprocess.run([compiler.compiler[0], "-dumpversion"], capture_output=True)
            if proc.returncode == 0:
                return semantic_version.Version(proc.stdout.decode())
    except:
        pass
    return semantic_version.Version("5.0.0") # assume lowest compatibility


def _set_cpp_flags(compiler, extension):
    # see `vendor/FAMSA/makefile`
    cc_version = _detect_compiler_version(compiler)

    # check we are not using the buggy clang version
    if _is_compiler_clang(compiler) and cc_version.major == 13 and cc_version.minor == 0:
        raise RuntimeError(
            "clang 13.0 is known to produce undefined behaviour "
            "(see https://github.com/refresh-bio/FAMSA/issues/32), "
            "please update or use a different compiler."
        )

    if cc_version.major <= 6:
        cpp_target = "11"
        defines = [("NO_PROFILE_PAR", 1), ("OLD_ATOMIC_FLAG", 1)]
    elif cc_version.major == 7:
        cpp_target = "14"
        defines = [("NO_PROFILE_PAR", 1), ("OLD_ATOMIC_FLAG", 1)]
    elif cc_version.major <= 10:
        cpp_target = "2a"
        defines = [("OLD_ATOMIC_FLAG", 1)]
    else:
        cpp_target = "20"
        defines = []

    extension.define_macros.extend(defines)
    if compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
        extension.extra_compile_args.append("-std=c++{}".format(cpp_target))
        extension.extra_compile_args.append("-pthread")
        extension.extra_compile_args.append("-fPIC")
        extension.extra_compile_args.append("-Wno-alloc-size-larger-than")
        extension.extra_compile_args.append("-Wno-char-subscripts")
        extension.extra_link_args.append("-pthread")
        extension.extra_link_args.append("-fPIC")
    elif compiler.compiler_type == "msvc":
        extension.extra_compile_args.append("/std:c{}".format(cpp_target))


def _set_debug_flags(compiler, extension):
    if compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
        extension.extra_compile_args.append("-g")
    elif compiler.compiler_type == "msvc":
        extension.extra_compile_args.append("/Z7")


# --- Extension with SIMD support --------------------------------------------

class Library(setuptools.extension.Library):

    def __init__(self, *args, **kwargs):
        self._needs_stub = False
        self.platform_sources = kwargs.pop("platform_sources", {})
        super().__init__(*args, **kwargs)


# --- Commands ------------------------------------------------------------------

class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    # --- Compatibility with `setuptools.Command`

    user_options = _build_ext.user_options + [
        ("disable-avx", None, "Force compiling the extension without AVX instructions"),
        ("disable-avx2", None, "Force compiling the extension without AVX2 instructions"),
        ("disable-neon", None, "Force compiling the extension without NEON instructions"),
    ]

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.disable_avx = False
        self.disable_avx2 = False
        self.disable_neon = False
        self.plat_name = sysconfig.get_platform()

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # detect platform options
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
        # force-disable AVX2 on MacOS (no OS support for runtime detection)
        if self.target_system == "macos":
            self.disable_avx2 = True
        # transfer arguments to the build_clib method
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.compiler = self.compiler
        self._clib_cmd.parallel = self.parallel
        self._clib_cmd.disable_avx = self.disable_avx
        self._clib_cmd.disable_avx2 = self.disable_avx2
        self._clib_cmd.disable_neon = self.disable_neon
        self._clib_cmd.plat_name = self.plat_name
        self._clib_cmd.target_machine = self.target_machine
        self._clib_cmd.target_system = self.target_system
        self._clib_cmd.target_cpu = self.target_cpu

    # --- Build code ---

    def build_extension(self, ext):
        # show the compiler being used
        _eprint("building", ext.name, "with", self.compiler.compiler_type, "compiler")

        # add debug symbols if we are building in debug mode
        if self.debug:
            _set_debug_flags(self.compiler, ext)
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))

        # add C++11 flags
        _set_cpp_flags(self.compiler, ext)

        # add Windows flags
        if self.target_system == "windows" and  self.compiler.compiler_type == "msvc":
            ext.define_macros.append(("WIN32", 1))

        # update link and include directories
        for name in ext.libraries:
            lib = self._clib_cmd.get_library(name)
            if lib is not None:
                libfile = self.compiler.library_filename(
                    lib.name, output_dir=self._clib_cmd.build_clib
                )
                ext.depends.append(libfile)
                ext.extra_objects.append(libfile)
                ext.extra_objects.extend(lib.extra_objects)
                ext.define_macros.extend(lib.define_macros)

        # build the rest of the extension as normal
        ext._needs_stub = False
        _build_ext.build_extension(self, ext)

    def build_extensions(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # remove universal compilation flags for OSX
        if self.target_system == "Darwin":
            _patch_osx_compiler(self.compiler, self.target_machine)

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "compiler_directives": {
                "cdivision": True,
                "nonecheck": False,
            },
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
                "DEFAULT_BUFFER_SIZE": io.DEFAULT_BUFFER_SIZE,
                "TARGET_CPU": self.target_cpu,
                "TARGET_SYSTEM": self.target_system,
            }
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["cdivision_warnings"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False

        # compile the C library
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # add the include dirs from the `build_clib` command
        for ext in self.extensions:
            for library in self.distribution.libraries:
                ext.include_dirs.extend(library.include_dirs)

        # cythonize the extensions and build normally
        self.extensions = cythonize(self.extensions, **cython_args)
        _build_ext.build_extensions(self)


class build_clib(_build_clib):
    """A custom `build_clib` that makes all C++ class attributes public.
    """

    # --- Compatibility with `setuptools.Command`

    user_options = _build_clib.user_options + [
        ("parallel=", "j", "number of parallel build jobs"),
        ("disable-avx", None, "Force compiling the library without AVX instructions"),
        ("disable-avx2", None, "Force compiling the library without AVX2 instructions"),
        ("disable-neon", None, "Force compiling the library without NEON instructions"),
    ]

    def initialize_options(self):
        _build_clib.initialize_options(self)
        self.parallel = None
        self.disable_avx = False
        self.disable_avx2 = False
        self.disable_neon = False
        self.plat_name = sysconfig.get_platform()

    def finalize_options(self):
        _build_clib.finalize_options(self)
        if self.parallel is not None:
            self.parallel = int(self.parallel)
        # detect platform options
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
        # force-disable AVX2 on MacOS (no OS support for runtime detection)
        if self.target_system == "macos":
            self.disable_avx2 = True
        # record SIMD-specific options
        self._simd_supported = dict(AVX=False, AVX2=False, NEON=False)
        self._simd_defines = dict(AVX=[], AVX2=[], NEON=[])
        self._simd_flags = dict(AVX=[], AVX2=[], NEON=[])
        self._simd_disabled = {
            "AVX": self.disable_avx,
            "AVX2": self.disable_avx2,
            "NEON": self.disable_neon,
        }

    # --- Autotools-like helpers ---

    def _patch_file(self, input, output):
        basename = os.path.basename(input)
        patchname = os.path.realpath(os.path.join(__file__, os.pardir, "patches", "{}.patch".format(basename)))
        _eprint("patching", os.path.relpath(input), "with", os.path.relpath(patchname))
        with open(patchname, "r") as patchfile:
            patch = patchfile.read()
        with open(input, "r") as src:
            srcdata = src.read()
        with open(output, "w") as dst:
            dst.write(_apply_patch(srcdata, patch))

    def _check_simd_generic(self, name, flags, header, vector, set, op, extract, result=1):
        _eprint('checking whether compiler can build', name, 'code', end="... ")

        base = "have_{}".format(name)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write("""
                #include <{}>
                int main() {{
                    {}      a = {}(1);
                            a = {}(a);
                    short   x = {}(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """.format(header, vector, set, op, extract))

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_preargs=flags)
            self.compiler.link_executable(objects, base, output_dir=self.build_temp)
            subprocess.run([binfile], check=True)
        except CompileError:
            _eprint("no")
            return False
        except (subprocess.SubprocessError, OSError):
            _eprint("yes, but cannot run code")
            return True  # assume we are cross-compiling, and still build
        else:
            if not flags:
                _eprint("yes")
            else:
                _eprint("yes, with {}".format(" ".join(flags)))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _avx_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:AVX2"]
        return ["-mavx", "-mpopcnt", "-funroll-loops"]

    def _check_avx(self):
        return self._check_simd_generic(
            "AVX",
            self._avx_flags(),
            header="immintrin.h",
            vector="__m256i",
            set="_mm256_set1_epi16",
            op="",
            extract="_mm256_extract_epi16",
        )

    def _avx2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:AVX2"]
        return ["-mavx2", "-mpopcnt", "-funroll-loops"]

    def _check_avx2(self):
        return self._check_simd_generic(
            "AVX2",
            self._avx2_flags(),
            header="immintrin.h",
            vector="__m256i",
            set="_mm256_set1_epi16",
            op="_mm256_abs_epi32",
            extract="_mm256_extract_epi16",
        )

    def _neon_flags(self):
        return ["-mfpu=neon"] if self.target_cpu == "arm" else []

    def _check_neon(self):
        return self._check_simd_generic(
            "NEON",
            self._neon_flags(),
            header="arm_neon.h",
            vector="int16x8_t",
            set="vdupq_n_s16",
            op="vabsq_s16",
            extract="vgetq_lane_s16"
        )

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    def get_library(self, name):
        return next((lib for lib in self.libraries if lib.name == name), None)

    # --- Build code ---

    def build_simd_code(self, library):
        # build platform-specific code
        for simd, sources in library.platform_sources.items():
            sources = [
                os.path.join(self.build_temp, "FAMSA", os.path.relpath(source, FAMSA_FOLDER))
                for source in sources
            ]
            if self._simd_supported[simd] and not self._simd_disabled[simd]:
                objects = [
                    s.replace(".cpp", self.compiler.obj_extension)
                    for s in sources
                ]
                for source, object in zip(sources, objects):
                    self.make_file(
                        [source],
                        object,
                        self.compiler.compile,
                        (
                            [source],
                            None,
                            library.define_macros + self._simd_defines[simd],
                            library.include_dirs,
                            self.debug,
                            library.extra_compile_args + self._simd_flags[simd],
                            None,
                            library.depends
                        )
                    )
                library.extra_objects.extend(objects)

    def build_libraries(self, libraries):
        # check `semantic_version` is available
        if isinstance(semantic_version, ImportError):
            raise RuntimeError("`semantic_version` is required to run `build_ext` command") from semantic_version

        # check for functions required for libcpu_features on OSX
        if self.target_system == "macos":
            _patch_osx_compiler(self.compiler, self.target_machine)

        # check if we can build platform-specific code
        if self.target_cpu == "x86":
            if not self._simd_disabled["AVX"] and self._check_avx():
                self._simd_supported["AVX"] = True
                self._simd_flags["AVX"].extend(self._avx_flags())
            if not self._simd_disabled["AVX2"] and self._check_avx2():
                self._simd_supported["AVX2"] = True
                self._simd_flags["AVX2"].extend(self._avx2_flags())
        elif self.target_cpu == "arm" or self.target_cpu == "aarch64":
            if not self._simd_disabled["NEON"] and self._check_neon():
                self._simd_supported["NEON"] = True
                self._simd_flags["NEON"].extend(self._neon_flags())

        # setup dispatcher of SIMD extensions
        for library in libraries:
            if self._simd_supported["AVX2"]:
                library.define_macros.append(("SIMD", 2))
            elif self._simd_supported["AVX"]:
                library.define_macros.append(("SIMD", 1))
            elif self._simd_supported["NEON"]:
                library.define_macros.append(("SIMD", 4))
            else:
                library.define_macros.append(("SIMD", 0))

        # build each library only if the sources are outdated
        self.mkpath(self.build_clib)
        for library in libraries:
            libname = self.compiler.library_filename(library.name, output_dir=self.build_clib)
            # self.make_file(library.sources, libname, self.build_library, (library,))
            self.build_library(library)

    def build_library(self, library):
        # show the compiler being used
        _eprint("building", library.name, "with", self.compiler.compiler_type, "compiler")

        # add debug flags if we are building in debug mode
        if self.debug:
            _set_debug_flags(self.compiler, library)
        # add C++11 flags
        if library.name == "famsa":
            _set_cpp_flags(self.compiler, library)

        # add Windows flags
        if self.target_system == "windows" and self.compiler.compiler_type == "msvc":
            library.define_macros.append(("WIN32", 1))

        # copy source code to build directory
        self.mkpath(os.path.join(self.build_clib, "FAMSA"))
        for dirpath, dirnames, filenames in os.walk(FAMSA_FOLDER):
            base = os.path.relpath(dirpath, FAMSA_FOLDER)
            self.mkpath(os.path.join(self.build_clib, "FAMSA", base))
            for filename in filenames:
                infile = os.path.join(dirpath, filename)
                outfile = os.path.join(self.build_clib, "FAMSA", base, filename)
                if os.path.exists(os.path.join(SETUP_FOLDER, "patches", "{}.patch".format(filename))):
                    self.make_file([infile], outfile, self._patch_file, (infile, outfile))
                else:
                    self.copy_file(infile, outfile)

        # fix library source paths
        for i, source in enumerate(library.sources):
            base = os.path.relpath(source, FAMSA_FOLDER)
            library.sources[i] = os.path.join(self.build_temp, "FAMSA", base)
        # fix include dirs to use the build directory folders
        for i, include_dir in enumerate(library.include_dirs):
            if include_dir.startswith(FAMSA_FOLDER):
                base = os.path.relpath(include_dir, FAMSA_FOLDER)
                library.include_dirs[i] = os.path.join(self.build_temp, "FAMSA", base)

        # store compile args
        compile_args = (
            None,
            library.define_macros,
            library.include_dirs + [self.build_clib],
            self.debug,
            library.extra_compile_args,
            None,
            library.depends,
        )
        # manually prepare sources and get the names of object files
        objects = [
            re.sub(r'(\.cpp|\.c)$', self.compiler.obj_extension, s)
            for s in library.sources
        ]
        # compile outdated files in parallel
        with multiprocessing.pool.ThreadPool(self.parallel) as pool:
            pool.starmap(
                functools.partial(self._compile_file, compile_args=compile_args),
                zip(library.sources, objects)
            )

        # build platform-specific code
        self.build_simd_code(library)

        # link into a static library
        libfile = self.compiler.library_filename(
            library.name,
            output_dir=self.build_clib,
        )
        self.make_file(
            objects,# + library.extra_objects,
            libfile,
            self.compiler.create_static_lib,
            (objects, library.name, self.build_clib, None, self.debug)
        )

    def _compile_file(self, source, object, compile_args):
        self.make_file(
            [source],
            object,
            self.compiler.compile,
            ([source], *compile_args)
        )


class clean(_clean):
    """A `clean` that removes intermediate files created by Cython.
    """

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pyfamsa")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c", "*.cpp"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                _eprint("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)


# --- Setup ---------------------------------------------------------------------

setuptools.setup(
    libraries=[
        # NOTE(@althonos): libdeflate is only needed for `IOService`, but we
        #                  don't use the FAMSA I/O so we don't have to build
        #                  it and link to it, as long as we're not building
        #                  `io_service.cpp` either
        Library(
            "famsa",
            language="c++",
            sources=[
                # COMMON_OBJS
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "msa.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "msa_refinement.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","AbstractTreeGenerator.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","Clustering.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","DistanceCalculator.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","FastTree.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","GuideTree.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","MSTPrim.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","NeighborJoining.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","NewickParser.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","SingleLinkage.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree","UPGMA.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "timer.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "log.cpp"),
                # os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "io_service.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "params.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "profile.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "profile_par.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "profile_seq.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "sequence.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "queues.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "libs", "mimalloc", "static.cpp"),
                # LCS_OBJS
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp.cpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_classic.cpp"),
            ],
            include_dirs=[
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "libs", "mimalloc"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "libs"),
                os.path.join(SETUP_FOLDER, "include"),
            ],
            depends=[
                # core
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "defs.h"),
                # os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "io_service.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "params.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "profile.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "queues.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "sequence.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "core", "version.h"),
                # lcs
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_avx_intr.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_avx2_intr.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_classic.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_neon_intr.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp.h"),
                # tree
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "AbstractTreeGenerator.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "AbstractTreeGenerator.hpp"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "Chained.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "Clustering.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "DistanceCalculator.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "FastTree.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "GuideTree.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "IPartialGenerator.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "MSTPrim.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "NeighborJoining.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "NewickParser.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "SingleLinkage.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "SingleLinkageQueue.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "TreeDefs.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "tree", "UPGMA.h"),
                # utils
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "array.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "conversion.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "cpuid.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "deterministic_random.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "log.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "memory_monotonic.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "meta_oper.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "pooled_threads.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "statistics.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "timer.h"),
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "utils.h"),
                # msa
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "msa.h"),
            ],
            platform_sources={
                "AVX": [
                    os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "utils_avx.cpp"),
                    os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_avx_intr.cpp"),
                ],
                "AVX2": [
                    os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "utils_avx2.cpp"),
                    os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_avx2_intr.cpp"),
                ],
                "NEON": [
                    os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "utils_neon.cpp"),
                    os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "lcs", "lcsbp_neon_intr.cpp"),
                ]
            },
        ),
    ],
    ext_modules=[
        Extension(
            "pyfamsa._famsa",
            language="c++",
            sources=[
                # for some reason, if this file is not compiled in the extension
                # this causes the `clear_mem` symbol to be missing from the
                # final shared object
                os.path.join(SETUP_FOLDER, "vendor", "FAMSA", "src", "utils", "utils.cpp"),
                os.path.join(SETUP_FOLDER, "pyfamsa", "_famsa.pyx"),
            ],
            include_dirs=[
                os.path.join(SETUP_FOLDER, "include"),
            ],
            libraries=[
                "famsa",
            ],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean
    }
)
