""" MSPtools """
from datetime import datetime

__version__ = '1.0.3'

def logger(message:str) -> None:
    """ Print message to the stdout. """
    print(
        f'[ {datetime.now().strftime(format="%Y-%m-%d %H:%M:%S")} ] {message}',
        flush=True
    )

def warning(msg:str) -> None:
    """ print a warning message """
    logger(f"[ MSPtools WARNING ] {msg}")
