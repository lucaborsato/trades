project('pytrades_lib', 'fortran',
        # default_options: ['warning_level=3', 'cpp_std=c++11', 'fortran_std=f2003'],
        # meson_version: '>=0.56.0'
        )

fortran_args = [
    '-cpp', 
    '-ffree-line-length-none', 
    '-Ofast', 
    '-fopenmp', 
    '-fPIC', 
    '-shared', 
]

# Aggiungi le opzioni di compilazione al compilatore Fortran
fortran_compiler = meson.get_compiler('fortran')
# fortran_compiler.set_options(fortran_args)

# Read package metadata from package.toml
# pkg_config = toml.load('package.toml')
# pkg_name = pkg_config['package']['name']
# pkg_version = pkg_config['package']['version']
# pkg_description = pkg_config['package']['description']
# pkg_author = pkg_config['package']['author']
# pkg_author_email = pkg_config['package']['author_email']
# pkg_license = pkg_config['package']['license']
# pkg_url = pkg_config['package']['url']
# pkg_classifiers = pkg_config['package']['classifiers']
# pkg_python_requires = pkg_config['package']['python_requires']

# Read Fortran source files from sources_list.txt
# fortran_sources = []
# sources_list = 'src/sources_list.txt'
# with open(sources_list) as f:
#     for line in f:
#         line = line.strip()
#         if line and not line.startswith('#'):
#             print(line)
#             fortran_sources += [line]
fortra_sources = [
    'src/constants.f90',
    'src/sort.f90',
    'src/custom_type.f90',
    'src/parameters.f90',
    'src/random_trades.f90',
    'src/convert_type.f90',
    'src/statistics.f90',
    'src/lin_fit.f90',
    'src/rotations.f90',
    'src/celestial_mechanics.f90',
    'src/eq_motion.f90',
    'src/parameters_conversion.f90',
    'src/linear_ephem.f90',
    'src/utils.f90',
    'src/init_trades.f90',
    'src/timing.f90',
    'src/gls_module.f90',
    'src/output_files.f90',
    'src/numerical_integrator.f90',
    'src/transits.f90',
    'src/radial_velocities.f90',
    'src/ode_run.f90',
    'src/derived_parameters_mod.f90',
    'src/fitness_module.f90',
    'src/grid_search.f90',
    'src/lm.f90',
    'src/pikaia.f90',
    'src/util_sort.f90',
    'src/util_qmc.f90',
    'src/opti_pso.f90',
    'src/gaussian.f90',
    'src/bootstrap.f90',
    'src/driver.f90'
]


# # Compile all Fortran files into a shared library
# fortran_lib = shared_library('pytrades_fortran',
#                              fortran_sources,
#                              include_directories: include_directories('src'),
#                              c_args: ['-cpp', '-ffree-line-length-none', '-Ofast', '-fopenmp', '-fPIC', '-shared', '-I.'],
#                              fortran_args: ['-cpp', '-ffree-line-length-none', '-Ofast', '-fopenmp', '-fPIC', '-shared', '-I.'],
#                              link_args: ['-lgomp'],
#                              install: true)
fortran_lib = shared_library('pytrades_fortran', fortra_sources,
    fortran_args: fortran_args,
    include_directories: include_directories('src')
)


# # Install the Python package
# subdir('pytrades')

# # pytrades/meson.build

# install_subdir('.', install_dir: join_paths(get_option('prefix'), 'lib', meson.project_name(), 'pytrades'))
