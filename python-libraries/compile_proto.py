import argparse
import os
import textwrap
from contextlib import contextmanager
from importlib import resources
from pathlib import Path


def handle_user_arguments(args=None) -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent(
        """\
    Compile the gRPC protocol files to python.
    """
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--proto-dir",
        required=True,
        dest="proto_dir",
        metavar="PATH",
        help="The path to the directory containing the proto files",
    )

    parser.add_argument(
        "--python-dir",
        required=True,
        dest="python_dir",
        metavar="PATH",
        help="The path to the directory where to generate the python files",
    )

    arguments = parser.parse_args(args)
    return arguments


def compile_protocol(proto_dir, python_dir):
    """
    Compile the protocol files to python.

    :param proto_dir: The path to the directory containing the proto files.
    :param python_dir: The path to the directory where to generate the python files.
    """
    from grpc_tools import protoc

    try:
        os.makedirs(python_dir)
    except OSError:
        pass

    # Note on calling grpc_tools.protoc as a python function:
    # grpc_tools.protoc.main is called by the command line with sys.argv and
    # the include path for the default protobuf proto files. sys.argv is a list of
    # the arguments passed to the command line, the first element of that list is
    # the command itself; what is passed as a command does not matter, but the actual
    # arguments must start at sys.argv[1] (hence "protoc" as first argument passed
    # to the function).
    proto_include = resources.files("grpc_tools") / "_proto"
    with move_in_directory(proto_dir):
        for protocol_file in Path(".").glob("**/*.proto"):
            print(f"Compiling {protocol_file}")
            protoc.main(
                (
                    "protoc",
                    "--proto_path=.",
                    f"--python_out={python_dir}",
                    f"--grpc_python_out={python_dir}",
                    str(protocol_file),
                    f"--proto_path={proto_include}",
                )
            )
    generated_protocol_directories = (
        path for path in (python_dir / "nanover/protocol").glob("**/*") if path.is_dir()
    )
    for directory in generated_protocol_directories:
        (directory / "__init__.py").touch()
        contained_files = (file for file in directory.glob("*_pb2*.py"))
        with open(directory / "__init__.py", "w+") as init_py:
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


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()

    compile_protocol(
        Path(arguments.proto_dir).resolve(),
        Path(arguments.python_dir).resolve(),
    )


if __name__ == "__main__":
    main()
