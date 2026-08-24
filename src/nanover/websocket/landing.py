from http.server import BaseHTTPRequestHandler, HTTPServer
from ssl import SSLContext

LANDING_HTML = """<dialog open style="inset: 0;"><button onclick="window.history.back(); window.close()" style="font-size: 10rem;">DONE</button></dialog>"""


class StaticHTTPRequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        self.wfile.write(LANDING_HTML.encode("UTF-8"))


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
