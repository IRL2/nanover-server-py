from http.server import HTTPServer, BaseHTTPRequestHandler
from ssl import SSLContext


class StaticHTTPRequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200, "HELLO")
        self.end_headers()
        self.wfile.write("HELLO".encode("UTF-8"))


def make_landing_page_server(
    *,
    host="0.0.0.0",
    port=0,
    ssl: SSLContext | None = None,
) -> HTTPServer:
    httpd = HTTPServer((host, port), StaticHTTPRequestHandler)

    if ssl:
        httpd.socket = ssl.wrap_socket(httpd.socket, server_side=True)

    return httpd
