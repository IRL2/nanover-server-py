from narupa.app import NarupaClient


def main():
    with NarupaClient() as client:
        frame = client.wait_until_first_frame()
        atoms = [index for index, id in enumerate(frame.arrays['atom.id']) if id == 'CG']
        print(atoms)


if __name__ == '__main__':
    main()
