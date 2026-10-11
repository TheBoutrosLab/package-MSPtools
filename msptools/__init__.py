""" MSPtools """
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version('msptools')
except PackageNotFoundError:
    __version__ = '0.1.0'

def logger(message:str) -> None:
    """ Print message to the stdout. """
    print(
        f'[ {datetime.now().strftime(format="%Y-%m-%d %H:%M:%S")} ] {message}',
        flush=True
    )

def warning(msg:str) -> None:
    """ print a warning message """
    logger(f"[ MSPtools WARNING ] {msg}")
