from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install
import sys

# NOTE: these's a lot of code here that's commented out, and that's
#       because it was the start of installing pyIMRPhenomD directly
#       from the MLDC repository. The module is, however, installed
#       if the MLDC master_install is done, so there's no need for
#       this as long as we stay in the docker. I'll leave the code
#       in case we change our mind though.

# these classes are required to pass additional arguments through pip
# to this script. Not entirely sure if everything in those classes is
# needed though, we might be able to shave off a few lines via trial
# and error.

#class CustomInstallCommand(install):
#    user_options = install.user_options + [('gsl=', None, None),
#                                           ('MLDC=', None, None)]
#    def initialize_options(self):
#        print('Running InstallCommand.initialize_options')
#        install.initialize_options(self)
#        self.gsl = None
#        self.MLDC = None
#
#    def finalize_options(self):
#        print('Running InstallCommand.finalize_options')
#        install.finalize_options(self)
#
#    def run(self):
#        print('Running InstallCommand.run')
#        install.run(self)
#
#
#class CustomDevelopCommand(develop):
#    user_options = develop.user_options + [('gsl=', None, None),
#                                           ('MLDC=', None, None)]
#    def initialize_options(self):
#        print('Running DevelopCommand.initialize_options')
#        develop.initialize_options(self)
#        self.gsl = None
#        self.MLDC = None
#
#    def finalize_options(self):
#        print('Running DevelopCommand.finalize_options')
#        develop.finalize_options(self)
#
#    def run(self):
#        print('Running DevelopCommand.run')
#        develop.run(self)

# the two extra commands needed for installing the IMRPhenomD
# waveform generator in the MLDC repository. We simply extract them
# from the argv list manually for now, and complain (without failing)
# if they are not found.

#gsl_location = None
#MLDC_location = None
## iterate over a copy, otherwise we would skip elements.
#for arg in sys.argv[:]:
#    if arg.startswith('--gsl='):
#        gsl_location = arg.split('=', 1)[1]
#        sys.argv.remove(arg)
#    if arg.startswith('--MLDC='):
#        MLDC_location = arg.split('=', 1)[1]
#        sys.argv.remove(arg)

# check if the required arguments were given, and skip making the waveform generator
# available if not. We may want to fail if these are not given eventually, but for
# now this seems more reasonable. We'd fail or default to the hdf5 route in the waveform
# generator interface.
#if gsl_location is None:
#    print('GSL location not given. IMR waveform generator cannot be installed.')
#if MLDC_location is None:
#    print('MLDC location not given. IMR waveform generator cannot be installed.')
#if gsl_location is not None and MLDC_location is not None:
#    raise NotImplementedError
#    # search / access IMRPD directory
#    # build .so
#    # copy .so into a local directory or write the required path to some file
#    # either import the file directly in the waveform generator interface, or add
#    #   the path and import from there


setup(name='jannas-utils',
      version='0.0',
      url='https://github.com/Janna112358/jannas-utils.git',
      author='Janna112358',
      packages=['jannas-utils'],
      zip_safe=False,
#      cmdclass={'install': CustomInstallCommand,
#                'develop': CustomDevelopCommand}
)
