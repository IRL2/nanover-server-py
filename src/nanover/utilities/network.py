def get_local_ip():
    import socket

    def attempt():
        yield from (
            ip
            for ip in socket.gethostbyname_ex(socket.gethostname())[2]
            if ip.startswith("192.")
        )

        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("8.8.8.8", 53))
            name = s.getsockname()[0]
            s.close()
        yield name

    ip = next(attempt())

    return ip
