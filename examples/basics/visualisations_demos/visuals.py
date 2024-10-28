from nanover.app.client import NanoverImdClient


def hiv_amprenavir():
    # Start a Narupa client, using connect_to_single_server with the ip address that started the server
    client = NanoverImdClient.autoconnect(name="my-nanover-server")
    client.subscribe_to_frames()

    # Wait for the first frame to be recieved about the current state of the system
    first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)

    print("Client connected, recieved frame with {0} atom(s)".format(first_frame.particle_count))

    my_first_selection = client.create_selection("Pro", list(range(0, 3129)))
    with my_first_selection.modify() as selection:
        selection.renderer = 'cartoon'

    my_second_selection = client.create_selection("Lig", list(range(3129, 3200)))

    with my_second_selection.modify() as selection:
        selection.renderer = {
            'color': {
                'type': 'cpk',
                'scheme': 'nanover'
            },
            'scale': 0.04,
            'render': 'liquorice'
        }
    print('Hiv cartoon, Amprenavir liquorice')
    client.close()
    return


def neuraminidase_zanamavir():
    # Start a Narupa client, using connect_to_single_server with the ip address that started the server
    client = NanoverImdClient.autoconnect(name="my-nanover-server")
    client.subscribe_to_frames()

    # Wait for the first frame to be recieved about the current state of the system
    first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)

    print("Client connected, recieved frame with {0} atom(s)".format(first_frame.particle_count))


    my_first_selection = client.create_selection("Pro", list(range(0, 5967)))
    with my_first_selection.modify() as selection:
        selection.renderer = 'cartoon'

    my_second_selection = client.create_selection("Lig", list(range(5967, 6009)))

    with my_second_selection.modify() as selection:
         selection.renderer = {
             'color': {
                 'type': 'cpk',
                 'scheme': 'nanover'
             },
             'scale': 0.04,
             'render': 'liquorice'
         }
    client.close()
    return


def trypsin_indoleamidine():
    # Start a Narupa client, using connect_to_single_server with the ip address that started the server
    client = NanoverImdClient.autoconnect(name="my-nanover-server")
    client.subscribe_to_frames()

    # Wait for the first frame to be recieved about the current state of the system
    first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)

    print("Client connected, recieved frame with {0} atom(s)".format(first_frame.particle_count))

    my_first_selection = client.create_selection("Pro", list(range(0, 3220)))
    with my_first_selection.modify() as selection:
        selection.renderer = 'cartoon'

    my_second_selection = client.create_selection("Lig", list(range(3220, 3257)))

    with my_second_selection.modify() as selection:
        selection.renderer = {
            'color': {
                'type': 'cpk',
                'scheme': 'nanover'
            },
            'scale': 0.04,
            'render': 'liquorice'
        }
    client.close()
    return

def clear_selections():
    # Start a Narupa client, using connect_to_single_server with the ip address that started the server
    client = NanoverImdClient.connect_to_single_server()
    client.subscribe_multiplayer()
    client.subscribe_to_frames()

    # Wait for the first frame to be recieved about the current state of the system
    first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)
    print("Client connected, recieved frame with {0} atom(s)".format(first_frame.particle_count))

    selection_ids = set(key for key in client._multiplayer_client.copy_state().keys() if key.startswith("selection."))
    for key in selection_ids:
        client.remove_shared_value(key)
    client.close()
    return
