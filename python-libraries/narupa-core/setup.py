# !/usr/bin/env python

# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.


import os
from contextlib import contextmanager
from pathlib import Path
import distutils.cmd
import distutils.log
from distutils.core import setup
from setuptools import find_namespace_packages


class CompileProtoCommand(distutils.cmd.Command):
    """
    Custom command to compile the protocol files.

    To run the command, call::

        python setup.py compile_proto

    While it is not recommended, it is possible to run the compilation with a
    different set of protocol files by passing the --proto-dir option:

        python setup.py compile_proto --proto-dir=../other-proto-files

    """

    description = "Compile the protocol files."
    user_options = [
        # The format is (long option, short option, description)
        ('proto-dir=', None, 'Path to the directory containing the protocol files.'),
    ]

    def initialize_options(self):
        """
        Set the default values for the options.
        """
        # By default, the setup.py file is in narupa-protocol/python-libraries/narupa-core,
        # the protocol files are in narupa-protocol/protocol, which is ../../protocol relative
        # to the directory of setup.py.
        here = Path(__file__).parent
        self.proto_dir = (here / '../../protocol')

    def finalize_options(self):
        """
        Post-process options.
        """
        self.proto_dir = Path(self.proto_dir)
        assert self.proto_dir.exists, 'The prototype directory {} does not exist.'.format(self.proto_dir)

    def run(self):
        """
        Run the compilation.
        """
        self.announce('Compile protocol files.', level=distutils.log.INFO)
        setup_path = Path(__file__).parent.resolve()
        compile_protocol(self.proto_dir, setup_path / 'src', self)


def compile_protocol(proto_dir, python_dir, logger):
    """
    Compile the protocol files to python.

    :param proto_dir: The path to the directory containing the proto files.
    :param python_dir: The path to the directory where to generate the python files.
    :param logger: The logger instance used by distutils.
    """
    from grpc_tools import protoc
    # Note on calling grpc_tools.protoc as a python function:
    # grpc_tools.protoc.main is called by the command line with sys.argv and
    # the include path for the default protobuf proto files. sys.argv is a list of
    # the arguments passed to the command line, the first element of that list is
    # the command itself; what is passed as a command does not matter, but the actual
    # arguments must start at sys.argv[1] (hence "protoc" as first argument passed
    # to the function).
    proto_include = protoc.pkg_resources.resource_filename('grpc_tools', '_proto')
    with move_in_directory(proto_dir):
        for protocol_file in Path('.').glob('**/*.proto'):
            logger.announce('Compiling {}'.format(protocol_file), level=distutils.log.INFO)
            protoc.main((
                'protoc',
                '--proto_path=.',
                '--python_out=' + str(python_dir),
                '--grpc_python_out=' + str(python_dir),
                str(protocol_file),
                '--proto_path={}'.format(proto_include),
            ))
    generated_protocol_directories = (path for path in (python_dir / 'narupa/protocol').glob('**/*') if path.is_dir())
    for directory in generated_protocol_directories:
        (directory / '__init__.py').touch()
        contained_files = (file for file in directory.glob('*_pb2*.py'))
        with open(directory / '__init__.py', "w+") as init_py:
            for contained_file in contained_files:
                file_name = os.path.splitext(os.path.split(contained_file)[1])[0]
                init_py.write("from .%s import *\n" % file_name)


@contextmanager
def move_in_directory(destination):
    """
    Context manager that moves in the given directory.

    When the interpreter enters the context manager, the working directory becomes
    the given destination. When the interpreter exists the context manager, the
    working directory is restored to where the working directory was before entering.

    Example:
    ========

    >>> with move_in_directory("bar"):
    >>>    # working directory is "bar"
    >>>    pass
    >>> # working directory is "foo" again

    :param destination: The directory to use as working directory.
    """
    destination = Path(destination)
    directory_to_restore = Path.cwd()
    try:
        os.chdir(str(destination))
        yield
    finally:
        os.chdir(str(directory_to_restore))


# Avoid repeating the content of requirements.txt here.
# The requirements.txt file is in the same directory as this setup.py,
# we then should now were the setup lies to be independent of the working
# directory.
requirements_path = Path(__file__).parent.resolve() / 'requirements.txt'
with open(str(requirements_path)) as f:
    requirements = f.readlines()

setup(name='narupa',
      version='1.0',
      description='Narupa python framework',
      author='Intangible Realities Lab',
      author_email='m.oconnor@bristol.ac.uk',
      url='https://gitlab.com/intangiblerealities/',
      packages=find_namespace_packages('src', include='narupa.*') + ['narupa.protocol'],
      package_dir={'': 'src'},
      package_data={
          '': ['py.typed']
      },
      install_requires=requirements,
      cmdclass={
          'compile_proto': CompileProtoCommand,
      },
      entry_points={
          'console_scripts': ['narupa-multiplayer=narupa.multiplayer.cli:main'],
      }
      )
