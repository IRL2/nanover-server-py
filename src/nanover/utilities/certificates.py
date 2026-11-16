import datetime
import ssl
from os import PathLike
from pathlib import Path

from cryptography import x509
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.x509.oid import NameOID
from platformdirs import PlatformDirs

DIRS = PlatformDirs("nanover-server", "irl")


def fetch_or_generate_self_signed_certificate_context(
    *,
    location: str | PathLike[str] | None = None,
    password=b"nanover",
):
    if location is None:
        location = DIRS.user_data_path
    else:
        location = Path(location)

    try:
        keyfile = location / "nanover-localhost-key.pem"
        certfile = location / "nanover-localhost-cert.pem"
        ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        ssl_context.load_verify_locations(certfile)
        ssl_context.load_cert_chain(certfile, keyfile=keyfile, password=password)
        return ssl_context
    except FileNotFoundError:
        generate_self_signed_certificate(location=location, password=password)
        return fetch_or_generate_self_signed_certificate_context(
            location=location, password=password
        )


def generate_self_signed_certificate(
    *,
    location: str | PathLike[str],
    password=b"nanover",
    duration=datetime.timedelta(weeks=52),
):
    """Generate a private key and self-signed certificate for localhost."""

    path = Path(location)
    path.mkdir(parents=True, exist_ok=True)

    key = rsa.generate_private_key(
        public_exponent=65537,
        key_size=2048,
    )

    with open(path / "nanover-localhost-key.pem", "wb") as f:
        f.write(
            key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.TraditionalOpenSSL,
                encryption_algorithm=serialization.BestAvailableEncryption(password),
            )
        )

    subject = issuer = x509.Name(
        [
            x509.NameAttribute(NameOID.COUNTRY_NAME, "XX"),
            x509.NameAttribute(NameOID.STATE_OR_PROVINCE_NAME, ""),
            x509.NameAttribute(NameOID.LOCALITY_NAME, ""),
            x509.NameAttribute(NameOID.ORGANIZATION_NAME, ""),
            x509.NameAttribute(NameOID.COMMON_NAME, "localhost"),
        ]
    )
    cert = (
        x509.CertificateBuilder()
        .subject_name(subject)
        .issuer_name(issuer)
        .public_key(key.public_key())
        .serial_number(x509.random_serial_number())
        .not_valid_before(datetime.datetime.now(datetime.UTC))
        .not_valid_after(datetime.datetime.now(datetime.UTC) + duration)
        .add_extension(
            x509.SubjectAlternativeName([x509.DNSName("localhost")]),
            critical=False,
        )
        .sign(key, hashes.SHA256())
    )
    with open(path / "nanover-localhost-cert.pem", "wb") as f:
        f.write(cert.public_bytes(serialization.Encoding.PEM))
