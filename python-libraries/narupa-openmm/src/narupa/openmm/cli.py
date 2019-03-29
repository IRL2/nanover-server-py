import argparse
from . import Runner


def handle_user_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--system', dest='system_xml_path')
    parser.add_argument('-s', '--structure', dest='structure_path')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    arguments = parser.parse_args()
    return arguments


def main():
    arguments = handle_user_arguments()
    simulation_runner = Runner.from_xml_input(
        input_xml=arguments.system_xml_path,
        pdb_path=arguments.structure_path,
    )
    simulation_runner.verbose = arguments.verbose
    simulation_runner.run()


if __name__ == '__main__':
    main()
