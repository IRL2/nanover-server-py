import argparse
import textwrap
from narupa.essd import DiscoveryClient, ServiceHub

LONG_TIME = 604800  # one week


def print_hub(hub: ServiceHub):
    services = ", ".join(f"{name} ({port})" for name, port in hub.services.items())
    print(f"{hub.name} ({hub.address}) -- {services}")


def handle_user_arguments(args=None) -> argparse.Namespace:
    description = textwrap.dedent("""\
    Run an ESSD client and print all services as they are found.
    """)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-p', '--port', default=None)
    parser.add_argument('-a', '--address', default=None)
    arguments = parser.parse_args(args)
    return arguments


def main():
    arguments = handle_user_arguments()

    try:
        with DiscoveryClient(arguments.address, arguments.port) as client:
            # search for LONG_TIME because there's no indefinite search
            # see: https://gitlab.com/intangiblerealities/narupa-protocol/issues/169
            for hub in client.search_for_services(LONG_TIME):
                print_hub(hub)
    except KeyboardInterrupt:
        print("Closing due to keyboard interrupt.")


if __name__ == '__main__':
    main()
